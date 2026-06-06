#define DUCKDB_EXTENSION_MAIN

#include "ducksmiles_extension.hpp"
#include "ducksmiles.h"

#include "duckdb.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "duckdb/common/vector_operations/unary_executor.hpp"
#include "duckdb/common/vector_operations/binary_executor.hpp"

#include <cmath>

namespace duckdb {

// ============================================================================
// Helper: call Rust string→bool/int FFI via UnaryExecutor
// ============================================================================

// VARCHAR → BOOLEAN (mol_is_valid, inchi_is_valid, inchikey_is_valid)
#define DEFINE_BOOL_FUNC(FuncName, RustFunc) \
static void FuncName(DataChunk &args, ExpressionState &state, Vector &result) { \
	UnaryExecutor::Execute<string_t, bool>(args.data[0], result, args.size(), \
		[](string_t input) -> bool { \
			return RustFunc((const uint8_t *)input.GetData(), input.GetSize()) == 1; \
		}); \
}

// VARCHAR → INTEGER with NULL on -1 sentinel
#define DEFINE_INT_FUNC(FuncName, RustFunc) \
static void FuncName(DataChunk &args, ExpressionState &state, Vector &result) { \
	UnaryExecutor::ExecuteWithNulls<string_t, int32_t>(args.data[0], result, args.size(), \
		[](string_t input, ValidityMask &mask, idx_t idx) -> int32_t { \
			int32_t val = RustFunc((const uint8_t *)input.GetData(), input.GetSize()); \
			if (val < 0) { mask.SetInvalid(idx); return 0; } \
			return val; \
		}); \
}

// VARCHAR → DOUBLE with NULL on NaN
#define DEFINE_DOUBLE_FUNC(FuncName, RustFunc) \
static void FuncName(DataChunk &args, ExpressionState &state, Vector &result) { \
	UnaryExecutor::ExecuteWithNulls<string_t, double>(args.data[0], result, args.size(), \
		[](string_t input, ValidityMask &mask, idx_t idx) -> double { \
			double val = RustFunc((const uint8_t *)input.GetData(), input.GetSize()); \
			if (std::isnan(val)) { mask.SetInvalid(idx); return 0.0; } \
			return val; \
		}); \
}

// VARCHAR → VARCHAR via Rust buffer-writing FFI. NULL on -1.
#define DEFINE_STR_FUNC(FuncName, RustFunc) \
static void FuncName(DataChunk &args, ExpressionState &state, Vector &result) { \
	UnaryExecutor::ExecuteWithNulls<string_t, string_t>(args.data[0], result, args.size(), \
		[&](string_t input, ValidityMask &mask, idx_t idx) -> string_t { \
			uint8_t buf[1024]; \
			int32_t len = RustFunc((const uint8_t *)input.GetData(), input.GetSize(), buf, sizeof(buf)); \
			if (len < 0) { mask.SetInvalid(idx); return string_t(); } \
			if (len == 0) { return StringVector::AddString(result, ""); } \
			return StringVector::AddString(result, (const char *)buf, len); \
		}); \
}

// ============================================================================
// SMILES functions
// ============================================================================

DEFINE_BOOL_FUNC(MolIsValidFunc, ds_mol_is_valid)
DEFINE_INT_FUNC(MolNumAtomsFunc, ds_mol_num_atoms)
DEFINE_INT_FUNC(MolNumBondsFunc, ds_mol_num_bonds)
DEFINE_STR_FUNC(MolFormulaFunc, ds_mol_formula)
DEFINE_DOUBLE_FUNC(MolWeightFunc, ds_mol_weight)
DEFINE_DOUBLE_FUNC(MolExactMassFunc, ds_mol_exact_mass)
DEFINE_DOUBLE_FUNC(LogpCrippenFunc, ds_logp_crippen)
DEFINE_DOUBLE_FUNC(TpsaFunc, ds_tpsa)

// add_hydrogens uses a larger 16KB buffer to handle drug-sized molecules
// (verbose SMILES with all H broken out can be ~5x the heavy-atom SMILES length).
static void AddHydrogensFunc(DataChunk &args, ExpressionState &state, Vector &result) {
	UnaryExecutor::ExecuteWithNulls<string_t, string_t>(args.data[0], result, args.size(),
		[&](string_t input, ValidityMask &mask, idx_t idx) -> string_t {
			uint8_t buf[16384];
			int32_t len = ds_add_hydrogens((const uint8_t *)input.GetData(), input.GetSize(), buf, sizeof(buf));
			if (len < 0) { mask.SetInvalid(idx); return string_t(); }
			if (len == 0) { return StringVector::AddString(result, ""); }
			return StringVector::AddString(result, (const char *)buf, len);
		});
}

// morgan_fp_bits(smi, radius, n_bits) → BLOB (ceil(n_bits/8) bytes).
// 16 KiB buffer is enough for up to 131072-bit fingerprints (standard sizes are
// 1024–4096).
static constexpr size_t MORGAN_BUF_BYTES = 16384;

static void MorganFpBitsFunc3(DataChunk &args, ExpressionState &state, Vector &result) {
	idx_t count = args.size();
	auto smi_data    = FlatVector::GetData<string_t>(args.data[0]);
	auto radius_data = FlatVector::GetData<int32_t>(args.data[1]);
	auto nbits_data  = FlatVector::GetData<int32_t>(args.data[2]);
	auto &validity   = FlatVector::Validity(result);
	for (idx_t i = 0; i < count; i++) {
		int32_t r = radius_data[i];
		int32_t n = nbits_data[i];
		if (r < 0 || n <= 0) { validity.SetInvalid(i); continue; }
		uint8_t buf[MORGAN_BUF_BYTES];
		int32_t len = ds_morgan_fp_bits(
			(const uint8_t *)smi_data[i].GetData(), smi_data[i].GetSize(),
			(uint32_t)r, (uint32_t)n, buf, sizeof(buf));
		if (len < 0) { validity.SetInvalid(i); continue; }
		auto blob = StringVector::AddStringOrBlob(result, (const char *)buf, len);
		FlatVector::GetData<string_t>(result)[i] = blob;
	}
}

// morgan_fp_bits(smi) → BLOB with defaults: ECFP4 (radius=2, 2048 bits → 256 bytes).
static void MorganFpBitsFunc1(DataChunk &args, ExpressionState &state, Vector &result) {
	UnaryExecutor::ExecuteWithNulls<string_t, string_t>(args.data[0], result, args.size(),
		[&](string_t input, ValidityMask &mask, idx_t idx) -> string_t {
			uint8_t buf[MORGAN_BUF_BYTES];
			int32_t len = ds_morgan_fp_bits(
				(const uint8_t *)input.GetData(), input.GetSize(),
				2u, 2048u, buf, sizeof(buf));
			if (len < 0) { mask.SetInvalid(idx); return string_t(); }
			return StringVector::AddStringOrBlob(result, (const char *)buf, len);
		});
}

// tanimoto_bit(BLOB, BLOB) → DOUBLE.
// Computes popcount(a & b) / popcount(a | b) on the raw BLOB bytes —
// no CAST AS BIT, no per-row table dispatch. Both arguments must be the
// same length (typically the 256-byte BLOB produced by morgan_fp_bits).
// BinaryExecutor handles flat / constant / dictionary vectors and NULL
// propagation; we only worry about length-mismatch (clear user error).
static void TanimotoBitFunc(DataChunk &args, ExpressionState &state, Vector &result) {
	BinaryExecutor::Execute<string_t, string_t, double>(
		args.data[0], args.data[1], result, args.size(),
		[](string_t a, string_t b) -> double {
			if (a.GetSize() != b.GetSize()) {
				throw InvalidInputException(
					"tanimoto_bit: BLOB lengths differ (%llu vs %llu bytes)",
					(unsigned long long)a.GetSize(),
					(unsigned long long)b.GetSize());
			}
			return ds_tanimoto_bit(
				(const uint8_t *)a.GetData(), a.GetSize(),
				(const uint8_t *)b.GetData(), b.GetSize());
		});
}

// ============================================================================
// InChI functions
// ============================================================================

DEFINE_BOOL_FUNC(InchiIsValidFunc, ds_inchi_is_valid)
DEFINE_BOOL_FUNC(InchiIsStandardFunc, ds_inchi_is_standard)
DEFINE_STR_FUNC(InchiVersionFunc, ds_inchi_version)
DEFINE_STR_FUNC(InchiFormulaFunc, ds_inchi_formula)
DEFINE_STR_FUNC(InchiConnectionsFunc, ds_inchi_connections)
DEFINE_STR_FUNC(InchiHydrogensFunc, ds_inchi_hydrogens)
DEFINE_STR_FUNC(InchiChargeFunc, ds_inchi_charge)
DEFINE_STR_FUNC(InchiStereoBondFunc, ds_inchi_stereo_bond)
DEFINE_STR_FUNC(InchiStereoTetrahedralFunc, ds_inchi_stereo_tetrahedral)
DEFINE_BOOL_FUNC(InchiHasStereoFunc, ds_inchi_has_stereo)
DEFINE_INT_FUNC(InchiNumStereoCentersFunc, ds_inchi_num_stereo_centers)

// ============================================================================
// InChIKey functions
// ============================================================================

DEFINE_BOOL_FUNC(InchikeyIsValidFunc, ds_inchikey_is_valid)
DEFINE_STR_FUNC(InchikeyConnectivityFunc, ds_inchikey_connectivity)
DEFINE_STR_FUNC(InchikeyStereofunc, ds_inchikey_stereo)
DEFINE_STR_FUNC(InchikeyProtonationFunc, ds_inchikey_protonation)

// inchi_skeleton_match(VARCHAR, VARCHAR) → BOOLEAN
static void InchiSkeletonMatchFunc(DataChunk &args, ExpressionState &state, Vector &result) {
	idx_t count = args.size();
	auto a_data = FlatVector::GetData<string_t>(args.data[0]);
	auto b_data = FlatVector::GetData<string_t>(args.data[1]);
	auto result_data = FlatVector::GetData<bool>(result);
	auto &validity = FlatVector::Validity(result);
	for (idx_t i = 0; i < count; i++) {
		int32_t val = ds_inchi_skeleton_match(
			(const uint8_t *)a_data[i].GetData(), a_data[i].GetSize(),
			(const uint8_t *)b_data[i].GetData(), b_data[i].GetSize());
		if (val < 0) {
			validity.SetInvalid(i);
			result_data[i] = false;
		} else {
			result_data[i] = (val == 1);
		}
	}
}

// ============================================================================
// MOL/SDF functions
// ============================================================================

DEFINE_STR_FUNC(MolBlockFormulaFunc, ds_mol_block_formula)
DEFINE_DOUBLE_FUNC(MolBlockWeightFunc, ds_mol_block_weight)
DEFINE_INT_FUNC(MolBlockNumAtomsFunc, ds_mol_block_num_atoms)
DEFINE_INT_FUNC(MolBlockNumBondsFunc, ds_mol_block_num_bonds)
DEFINE_STR_FUNC(MolBlockNameFunc, ds_mol_block_name)
DEFINE_INT_FUNC(SdfCountFunc, ds_sdf_count)

// ============================================================================
// PDB/CIF/XYZ functions — format arg hardcoded to 0 (auto-detect)
// ============================================================================

#define DEFINE_STRUCTURE_INT_FUNC(FuncName, RustFunc) \
static void FuncName(DataChunk &args, ExpressionState &state, Vector &result) { \
	UnaryExecutor::ExecuteWithNulls<string_t, int32_t>(args.data[0], result, args.size(), \
		[](string_t input, ValidityMask &mask, idx_t idx) -> int32_t { \
			int32_t val = RustFunc((const uint8_t *)input.GetData(), input.GetSize(), 0); \
			if (val < 0) { mask.SetInvalid(idx); return 0; } \
			return val; \
		}); \
}

DEFINE_STRUCTURE_INT_FUNC(StructureAtomCountFunc, ds_structure_atom_count)
DEFINE_STRUCTURE_INT_FUNC(StructureChainCountFunc, ds_structure_chain_count)
DEFINE_STRUCTURE_INT_FUNC(StructureResidueCountFunc, ds_structure_residue_count)
DEFINE_STRUCTURE_INT_FUNC(StructureModelCountFunc, ds_structure_model_count)

// ============================================================================
// SELFIES functions
// ============================================================================

DEFINE_STR_FUNC(SmilesToSelfiesFunc, ds_smiles_to_selfies)
DEFINE_STR_FUNC(SelfiesToSmilesFunc, ds_selfies_to_smiles)
DEFINE_BOOL_FUNC(SelfiesIsValidFunc, ds_selfies_is_valid)

// ============================================================================
// Registration
// ============================================================================

static void RegisterDucksmilesFunctions(ExtensionLoader &loader) {
	// --- SMILES ---
	loader.RegisterFunction(ScalarFunction("mol_is_valid",    {LogicalType::VARCHAR}, LogicalType::BOOLEAN, MolIsValidFunc));
	loader.RegisterFunction(ScalarFunction("mol_formula",     {LogicalType::VARCHAR}, LogicalType::VARCHAR, MolFormulaFunc));
	loader.RegisterFunction(ScalarFunction("mol_num_atoms",   {LogicalType::VARCHAR}, LogicalType::INTEGER, MolNumAtomsFunc));
	loader.RegisterFunction(ScalarFunction("mol_num_bonds",   {LogicalType::VARCHAR}, LogicalType::INTEGER, MolNumBondsFunc));
	loader.RegisterFunction(ScalarFunction("mol_weight",      {LogicalType::VARCHAR}, LogicalType::DOUBLE,  MolWeightFunc));
	loader.RegisterFunction(ScalarFunction("mol_exact_mass",  {LogicalType::VARCHAR}, LogicalType::DOUBLE,  MolExactMassFunc));
	loader.RegisterFunction(ScalarFunction("logp_crippen",    {LogicalType::VARCHAR}, LogicalType::DOUBLE,  LogpCrippenFunc));
	loader.RegisterFunction(ScalarFunction("tpsa",            {LogicalType::VARCHAR}, LogicalType::DOUBLE,  TpsaFunc));
	loader.RegisterFunction(ScalarFunction("add_hydrogens",   {LogicalType::VARCHAR}, LogicalType::VARCHAR, AddHydrogensFunc));
	loader.RegisterFunction(ScalarFunction("morgan_fp_bits",  {LogicalType::VARCHAR}, LogicalType::BLOB,    MorganFpBitsFunc1));
	loader.RegisterFunction(ScalarFunction("morgan_fp_bits",
		{LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::INTEGER},
		LogicalType::BLOB, MorganFpBitsFunc3));
	loader.RegisterFunction(ScalarFunction("tanimoto_bit",
		{LogicalType::BLOB, LogicalType::BLOB},
		LogicalType::DOUBLE, TanimotoBitFunc));

	// --- InChI layer extraction ---
	loader.RegisterFunction(ScalarFunction("inchi_is_valid",           {LogicalType::VARCHAR}, LogicalType::BOOLEAN, InchiIsValidFunc));
	loader.RegisterFunction(ScalarFunction("inchi_is_standard",        {LogicalType::VARCHAR}, LogicalType::BOOLEAN, InchiIsStandardFunc));
	loader.RegisterFunction(ScalarFunction("inchi_version",            {LogicalType::VARCHAR}, LogicalType::VARCHAR, InchiVersionFunc));
	loader.RegisterFunction(ScalarFunction("inchi_formula",            {LogicalType::VARCHAR}, LogicalType::VARCHAR, InchiFormulaFunc));
	loader.RegisterFunction(ScalarFunction("inchi_connections",        {LogicalType::VARCHAR}, LogicalType::VARCHAR, InchiConnectionsFunc));
	loader.RegisterFunction(ScalarFunction("inchi_hydrogens",          {LogicalType::VARCHAR}, LogicalType::VARCHAR, InchiHydrogensFunc));
	loader.RegisterFunction(ScalarFunction("inchi_charge",             {LogicalType::VARCHAR}, LogicalType::VARCHAR, InchiChargeFunc));
	loader.RegisterFunction(ScalarFunction("inchi_stereo_bond",        {LogicalType::VARCHAR}, LogicalType::VARCHAR, InchiStereoBondFunc));
	loader.RegisterFunction(ScalarFunction("inchi_stereo_tetrahedral", {LogicalType::VARCHAR}, LogicalType::VARCHAR, InchiStereoTetrahedralFunc));
	loader.RegisterFunction(ScalarFunction("inchi_has_stereo",         {LogicalType::VARCHAR}, LogicalType::BOOLEAN, InchiHasStereoFunc));
	loader.RegisterFunction(ScalarFunction("inchi_num_stereo_centers", {LogicalType::VARCHAR}, LogicalType::INTEGER, InchiNumStereoCentersFunc));
	// --- InChIKey ---
	loader.RegisterFunction(ScalarFunction("inchikey_is_valid",      {LogicalType::VARCHAR}, LogicalType::BOOLEAN, InchikeyIsValidFunc));
	loader.RegisterFunction(ScalarFunction("inchikey_connectivity",  {LogicalType::VARCHAR}, LogicalType::VARCHAR, InchikeyConnectivityFunc));
	loader.RegisterFunction(ScalarFunction("inchikey_stereo",        {LogicalType::VARCHAR}, LogicalType::VARCHAR, InchikeyStereofunc));
	loader.RegisterFunction(ScalarFunction("inchikey_protonation",   {LogicalType::VARCHAR}, LogicalType::VARCHAR, InchikeyProtonationFunc));

	// --- Comparison ---
	loader.RegisterFunction(ScalarFunction("inchi_skeleton_match",   {LogicalType::VARCHAR, LogicalType::VARCHAR}, LogicalType::BOOLEAN, InchiSkeletonMatchFunc));

	// --- MOL/SDF ---
	loader.RegisterFunction(ScalarFunction("mol_block_formula",    {LogicalType::VARCHAR}, LogicalType::VARCHAR, MolBlockFormulaFunc));
	loader.RegisterFunction(ScalarFunction("mol_block_weight",     {LogicalType::VARCHAR}, LogicalType::DOUBLE,  MolBlockWeightFunc));
	loader.RegisterFunction(ScalarFunction("mol_block_num_atoms",  {LogicalType::VARCHAR}, LogicalType::INTEGER, MolBlockNumAtomsFunc));
	loader.RegisterFunction(ScalarFunction("mol_block_num_bonds",  {LogicalType::VARCHAR}, LogicalType::INTEGER, MolBlockNumBondsFunc));
	loader.RegisterFunction(ScalarFunction("mol_block_name",       {LogicalType::VARCHAR}, LogicalType::VARCHAR, MolBlockNameFunc));
	loader.RegisterFunction(ScalarFunction("sdf_count",            {LogicalType::VARCHAR}, LogicalType::INTEGER, SdfCountFunc));

	// --- PDB/CIF/XYZ structure ---
	loader.RegisterFunction(ScalarFunction("structure_atom_count",    {LogicalType::VARCHAR}, LogicalType::INTEGER, StructureAtomCountFunc));
	loader.RegisterFunction(ScalarFunction("structure_chain_count",   {LogicalType::VARCHAR}, LogicalType::INTEGER, StructureChainCountFunc));
	loader.RegisterFunction(ScalarFunction("structure_residue_count", {LogicalType::VARCHAR}, LogicalType::INTEGER, StructureResidueCountFunc));
	loader.RegisterFunction(ScalarFunction("structure_model_count",   {LogicalType::VARCHAR}, LogicalType::INTEGER, StructureModelCountFunc));

	// --- SELFIES ---
	loader.RegisterFunction(ScalarFunction("smiles_to_selfies",    {LogicalType::VARCHAR}, LogicalType::VARCHAR, SmilesToSelfiesFunc));
	loader.RegisterFunction(ScalarFunction("selfies_to_smiles",    {LogicalType::VARCHAR}, LogicalType::VARCHAR, SelfiesToSmilesFunc));
	loader.RegisterFunction(ScalarFunction("selfies_is_valid",     {LogicalType::VARCHAR}, LogicalType::BOOLEAN, SelfiesIsValidFunc));
}

// ============================================================================
// Extension lifecycle
// ============================================================================

void DucksmilesExtension::Load(ExtensionLoader &loader) {
	RegisterDucksmilesFunctions(loader);
}

std::string DucksmilesExtension::Name() {
	return "ducksmiles";
}

std::string DucksmilesExtension::Version() const {
#ifdef EXT_VERSION_DUCKSMILES
	return EXT_VERSION_DUCKSMILES;
#else
	return "0.1.0";
#endif
}

} // namespace duckdb

extern "C" {

DUCKDB_EXTENSION_API void ducksmiles_duckdb_cpp_init(duckdb::ExtensionLoader &loader) {
	duckdb::DucksmilesExtension ext;
	ext.Load(loader);
}

DUCKDB_EXTENSION_API void ducksmiles_init(duckdb::DatabaseInstance &db) {
	// Legacy init
}

DUCKDB_EXTENSION_API const char *ducksmiles_version() {
	return duckdb::DuckDB::LibraryVersion();
}

}
