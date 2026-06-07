#define DUCKDB_EXTENSION_MAIN

#include "ducksmiles_extension.hpp"
#include "ducksmiles.h"

#include "duckdb.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "duckdb/common/vector_operations/unary_executor.hpp"
#include "duckdb/common/vector_operations/binary_executor.hpp"

#include <cmath>
#include <vector>

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

// VARCHAR → BOOLEAN with NULL on -1 sentinel
#define DEFINE_SENTINEL_BOOL_FUNC(FuncName, RustFunc) \
static void FuncName(DataChunk &args, ExpressionState &state, Vector &result) { \
	UnaryExecutor::ExecuteWithNulls<string_t, bool>(args.data[0], result, args.size(), \
		[](string_t input, ValidityMask &mask, idx_t idx) -> bool { \
			int32_t val = RustFunc((const uint8_t *)input.GetData(), input.GetSize()); \
			if (val < 0) { mask.SetInvalid(idx); return false; } \
			return val == 1; \
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

template <class Callback>
static string_t DynamicStringResult(Vector &result, ValidityMask &mask, idx_t idx, Callback callback) {
	int32_t needed = callback(nullptr, 0);
	if (needed < 0) {
		mask.SetInvalid(idx);
		return string_t();
	}
	if (needed == 0) {
		return StringVector::AddString(result, "");
	}

	std::vector<uint8_t> buf((size_t)needed);
	int32_t actual = callback(buf.data(), buf.size());
	if (actual < 0) {
		mask.SetInvalid(idx);
		return string_t();
	}
	if (actual != needed) {
		throw InvalidInputException("ducksmiles string FFI length changed between sizing and writing");
	}
	return StringVector::AddString(result, (const char *)buf.data(), buf.size());
}

#define DEFINE_DYNAMIC_STR_FUNC(FuncName, RustFunc) \
static void FuncName(DataChunk &args, ExpressionState &state, Vector &result) { \
	UnaryExecutor::ExecuteWithNulls<string_t, string_t>(args.data[0], result, args.size(), \
		[&](string_t input, ValidityMask &mask, idx_t idx) -> string_t { \
			return DynamicStringResult(result, mask, idx, [&](uint8_t *out, size_t cap) -> int32_t { \
				return RustFunc((const uint8_t *)input.GetData(), input.GetSize(), out, cap); \
			}); \
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
DEFINE_STR_FUNC(CanonicalSmilesFunc, ds_canonical_smiles)
DEFINE_DYNAMIC_STR_FUNC(MurckoScaffoldFunc, ds_murcko_scaffold)
DEFINE_DYNAMIC_STR_FUNC(GenericScaffoldFunc, ds_generic_scaffold)
DEFINE_DYNAMIC_STR_FUNC(RingSystemsJsonFunc, ds_ring_systems_json)
DEFINE_INT_FUNC(NumHAcceptorsFunc, ds_num_h_acceptors)
DEFINE_INT_FUNC(NumHDonorsFunc, ds_num_h_donors)
DEFINE_INT_FUNC(NumRotatableBondsFunc, ds_num_rotatable_bonds)
DEFINE_INT_FUNC(RingCountFunc, ds_ring_count)
DEFINE_INT_FUNC(NumAromaticRingsFunc, ds_num_aromatic_rings)
DEFINE_INT_FUNC(NumHeteroatomsFunc, ds_num_heteroatoms)
DEFINE_DOUBLE_FUNC(FractionCsp3Func, ds_fraction_csp3)

static void MolHasSubstructureFunc(DataChunk &args, ExpressionState &state, Vector &result) {
	idx_t count = args.size();
	auto smi_data = FlatVector::GetData<string_t>(args.data[0]);
	auto smarts_data = FlatVector::GetData<string_t>(args.data[1]);
	auto result_data = FlatVector::GetData<bool>(result);
	auto &validity = FlatVector::Validity(result);
	for (idx_t i = 0; i < count; i++) {
		int32_t val = ds_mol_has_substructure(
			(const uint8_t *)smi_data[i].GetData(), smi_data[i].GetSize(),
			(const uint8_t *)smarts_data[i].GetData(), smarts_data[i].GetSize());
		if (val < 0) {
			validity.SetInvalid(i);
			result_data[i] = false;
		} else {
			result_data[i] = val == 1;
		}
	}
}

static void MolSubstructureCountFunc(DataChunk &args, ExpressionState &state, Vector &result) {
	idx_t count = args.size();
	auto smi_data = FlatVector::GetData<string_t>(args.data[0]);
	auto smarts_data = FlatVector::GetData<string_t>(args.data[1]);
	auto result_data = FlatVector::GetData<int32_t>(result);
	auto &validity = FlatVector::Validity(result);
	for (idx_t i = 0; i < count; i++) {
		int32_t val = ds_mol_substructure_count(
			(const uint8_t *)smi_data[i].GetData(), smi_data[i].GetSize(),
			(const uint8_t *)smarts_data[i].GetData(), smarts_data[i].GetSize());
		if (val < 0) {
			validity.SetInvalid(i);
			result_data[i] = 0;
		} else {
			result_data[i] = val;
		}
	}
}

static void MolSubstructureMatchesJsonFunc(DataChunk &args, ExpressionState &state, Vector &result) {
	idx_t count = args.size();
	args.data[0].Flatten(count);
	args.data[1].Flatten(count);
	auto smi_data = FlatVector::GetData<string_t>(args.data[0]);
	auto smarts_data = FlatVector::GetData<string_t>(args.data[1]);
	auto result_data = FlatVector::GetData<string_t>(result);
	auto &validity = FlatVector::Validity(result);
	auto &smi_validity = FlatVector::Validity(args.data[0]);
	auto &smarts_validity = FlatVector::Validity(args.data[1]);
	for (idx_t i = 0; i < count; i++) {
		if (!smi_validity.RowIsValid(i) || !smarts_validity.RowIsValid(i)) {
			validity.SetInvalid(i);
			continue;
		}
		auto smi = smi_data[i];
		auto smarts = smarts_data[i];
		result_data[i] = DynamicStringResult(result, validity, i, [&](uint8_t *out, size_t cap) -> int32_t {
			return ds_mol_substructure_matches_json(
				(const uint8_t *)smi.GetData(), smi.GetSize(),
				(const uint8_t *)smarts.GetData(), smarts.GetSize(),
				out, cap);
		});
	}
}

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

// maccs_keys(smi) → BLOB (fixed 21 bytes, 167-bit MACCS vector).
static constexpr size_t MACCS_BUF_BYTES = 21;

static void MaccsKeysFunc(DataChunk &args, ExpressionState &state, Vector &result) {
	UnaryExecutor::ExecuteWithNulls<string_t, string_t>(args.data[0], result, args.size(),
		[&](string_t input, ValidityMask &mask, idx_t idx) -> string_t {
			uint8_t buf[MACCS_BUF_BYTES];
			int32_t len = ds_maccs_keys(
				(const uint8_t *)input.GetData(), input.GetSize(),
				buf, sizeof(buf));
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
DEFINE_DYNAMIC_STR_FUNC(MolBlockPropertiesJsonFunc, ds_mol_block_properties_json)
DEFINE_DYNAMIC_STR_FUNC(MolBlockAtomsJsonFunc, ds_mol_block_atoms_json)
DEFINE_DYNAMIC_STR_FUNC(MolBlockBondsJsonFunc, ds_mol_block_bonds_json)
DEFINE_DYNAMIC_STR_FUNC(MolBlockJsonFunc, ds_mol_block_json)
DEFINE_DYNAMIC_STR_FUNC(SdfPropertiesJsonFunc, ds_sdf_properties_json)
DEFINE_INT_FUNC(SdfCountFunc, ds_sdf_count)
DEFINE_SENTINEL_BOOL_FUNC(MolBlockHas3dFunc, ds_mol_block_has_3d)
DEFINE_DOUBLE_FUNC(MolBlockCentroidXFunc, ds_mol_block_centroid_x)
DEFINE_DOUBLE_FUNC(MolBlockCentroidYFunc, ds_mol_block_centroid_y)
DEFINE_DOUBLE_FUNC(MolBlockCentroidZFunc, ds_mol_block_centroid_z)
DEFINE_DOUBLE_FUNC(MolBlockRadiusOfGyrationFunc, ds_mol_block_radius_of_gyration)
DEFINE_DOUBLE_FUNC(MolBlockMinXFunc, ds_mol_block_min_x)
DEFINE_DOUBLE_FUNC(MolBlockMaxXFunc, ds_mol_block_max_x)
DEFINE_DOUBLE_FUNC(MolBlockMinYFunc, ds_mol_block_min_y)
DEFINE_DOUBLE_FUNC(MolBlockMaxYFunc, ds_mol_block_max_y)
DEFINE_DOUBLE_FUNC(MolBlockMinZFunc, ds_mol_block_min_z)
DEFINE_DOUBLE_FUNC(MolBlockMaxZFunc, ds_mol_block_max_z)

static void MolBlockPropertyFunc(DataChunk &args, ExpressionState &state, Vector &result) {
	idx_t count = args.size();
	args.data[0].Flatten(count);
	args.data[1].Flatten(count);
	auto mol_data = FlatVector::GetData<string_t>(args.data[0]);
	auto key_data = FlatVector::GetData<string_t>(args.data[1]);
	auto result_data = FlatVector::GetData<string_t>(result);
	auto &validity = FlatVector::Validity(result);
	auto &mol_validity = FlatVector::Validity(args.data[0]);
	auto &key_validity = FlatVector::Validity(args.data[1]);
	for (idx_t i = 0; i < count; i++) {
		if (!mol_validity.RowIsValid(i) || !key_validity.RowIsValid(i)) {
			validity.SetInvalid(i);
			continue;
		}
		auto mol = mol_data[i];
		auto key = key_data[i];
		result_data[i] = DynamicStringResult(result, validity, i, [&](uint8_t *out, size_t cap) -> int32_t {
			return ds_mol_block_property(
				(const uint8_t *)mol.GetData(), mol.GetSize(),
				(const uint8_t *)key.GetData(), key.GetSize(),
				out, cap);
		});
	}
}

static void SdfPropertyFunc(DataChunk &args, ExpressionState &state, Vector &result) {
	idx_t count = args.size();
	args.data[0].Flatten(count);
	args.data[1].Flatten(count);
	args.data[2].Flatten(count);
	auto sdf_data = FlatVector::GetData<string_t>(args.data[0]);
	auto index_data = FlatVector::GetData<int32_t>(args.data[1]);
	auto key_data = FlatVector::GetData<string_t>(args.data[2]);
	auto result_data = FlatVector::GetData<string_t>(result);
	auto &validity = FlatVector::Validity(result);
	auto &sdf_validity = FlatVector::Validity(args.data[0]);
	auto &index_validity = FlatVector::Validity(args.data[1]);
	auto &key_validity = FlatVector::Validity(args.data[2]);
	for (idx_t i = 0; i < count; i++) {
		if (!sdf_validity.RowIsValid(i) || !index_validity.RowIsValid(i) || !key_validity.RowIsValid(i)) {
			validity.SetInvalid(i);
			continue;
		}
		auto sdf = sdf_data[i];
		auto key = key_data[i];
		result_data[i] = DynamicStringResult(result, validity, i, [&](uint8_t *out, size_t cap) -> int32_t {
			return ds_sdf_property(
				(const uint8_t *)sdf.GetData(), sdf.GetSize(),
				index_data[i],
				(const uint8_t *)key.GetData(), key.GetSize(),
				out, cap);
		});
	}
}

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

#define DEFINE_STRUCTURE_DOUBLE_FUNC(FuncName, RustFunc) \
static void FuncName(DataChunk &args, ExpressionState &state, Vector &result) { \
	UnaryExecutor::ExecuteWithNulls<string_t, double>(args.data[0], result, args.size(), \
		[](string_t input, ValidityMask &mask, idx_t idx) -> double { \
			double val = RustFunc((const uint8_t *)input.GetData(), input.GetSize(), 0); \
			if (std::isnan(val)) { mask.SetInvalid(idx); return 0.0; } \
			return val; \
		}); \
}

DEFINE_STRUCTURE_DOUBLE_FUNC(StructureCentroidXFunc, ds_structure_centroid_x)
DEFINE_STRUCTURE_DOUBLE_FUNC(StructureCentroidYFunc, ds_structure_centroid_y)
DEFINE_STRUCTURE_DOUBLE_FUNC(StructureCentroidZFunc, ds_structure_centroid_z)
DEFINE_STRUCTURE_DOUBLE_FUNC(StructureRadiusOfGyrationFunc, ds_structure_radius_of_gyration)
DEFINE_STRUCTURE_DOUBLE_FUNC(StructureMinXFunc, ds_structure_min_x)
DEFINE_STRUCTURE_DOUBLE_FUNC(StructureMaxXFunc, ds_structure_max_x)
DEFINE_STRUCTURE_DOUBLE_FUNC(StructureMinYFunc, ds_structure_min_y)
DEFINE_STRUCTURE_DOUBLE_FUNC(StructureMaxYFunc, ds_structure_max_y)
DEFINE_STRUCTURE_DOUBLE_FUNC(StructureMinZFunc, ds_structure_min_z)
DEFINE_STRUCTURE_DOUBLE_FUNC(StructureMaxZFunc, ds_structure_max_z)

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
	loader.RegisterFunction(ScalarFunction("canonical_smiles", {LogicalType::VARCHAR}, LogicalType::VARCHAR, CanonicalSmilesFunc));
	loader.RegisterFunction(ScalarFunction("murcko_scaffold", {LogicalType::VARCHAR}, LogicalType::VARCHAR, MurckoScaffoldFunc));
	loader.RegisterFunction(ScalarFunction("generic_scaffold", {LogicalType::VARCHAR}, LogicalType::VARCHAR, GenericScaffoldFunc));
	loader.RegisterFunction(ScalarFunction("ring_systems_json", {LogicalType::VARCHAR}, LogicalType::VARCHAR, RingSystemsJsonFunc));
	loader.RegisterFunction(ScalarFunction("num_h_acceptors", {LogicalType::VARCHAR}, LogicalType::INTEGER, NumHAcceptorsFunc));
	loader.RegisterFunction(ScalarFunction("num_h_donors",    {LogicalType::VARCHAR}, LogicalType::INTEGER, NumHDonorsFunc));
	loader.RegisterFunction(ScalarFunction("num_rotatable_bonds", {LogicalType::VARCHAR}, LogicalType::INTEGER, NumRotatableBondsFunc));
	loader.RegisterFunction(ScalarFunction("ring_count",      {LogicalType::VARCHAR}, LogicalType::INTEGER, RingCountFunc));
	loader.RegisterFunction(ScalarFunction("num_aromatic_rings", {LogicalType::VARCHAR}, LogicalType::INTEGER, NumAromaticRingsFunc));
	loader.RegisterFunction(ScalarFunction("num_heteroatoms", {LogicalType::VARCHAR}, LogicalType::INTEGER, NumHeteroatomsFunc));
	loader.RegisterFunction(ScalarFunction("fraction_csp3",   {LogicalType::VARCHAR}, LogicalType::DOUBLE,  FractionCsp3Func));
	loader.RegisterFunction(ScalarFunction("mol_has_substructure",
		{LogicalType::VARCHAR, LogicalType::VARCHAR},
		LogicalType::BOOLEAN, MolHasSubstructureFunc));
	loader.RegisterFunction(ScalarFunction("mol_substructure_count",
		{LogicalType::VARCHAR, LogicalType::VARCHAR},
		LogicalType::INTEGER, MolSubstructureCountFunc));
	loader.RegisterFunction(ScalarFunction("mol_substructure_matches_json",
		{LogicalType::VARCHAR, LogicalType::VARCHAR},
		LogicalType::VARCHAR, MolSubstructureMatchesJsonFunc));
	loader.RegisterFunction(ScalarFunction("add_hydrogens",   {LogicalType::VARCHAR}, LogicalType::VARCHAR, AddHydrogensFunc));
	loader.RegisterFunction(ScalarFunction("morgan_fp_bits",  {LogicalType::VARCHAR}, LogicalType::BLOB,    MorganFpBitsFunc1));
	loader.RegisterFunction(ScalarFunction("morgan_fp_bits",
		{LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::INTEGER},
		LogicalType::BLOB, MorganFpBitsFunc3));
	loader.RegisterFunction(ScalarFunction("maccs_keys",      {LogicalType::VARCHAR}, LogicalType::BLOB,    MaccsKeysFunc));
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
	loader.RegisterFunction(ScalarFunction("mol_block_property",   {LogicalType::VARCHAR, LogicalType::VARCHAR}, LogicalType::VARCHAR, MolBlockPropertyFunc));
	loader.RegisterFunction(ScalarFunction("mol_block_properties_json", {LogicalType::VARCHAR}, LogicalType::VARCHAR, MolBlockPropertiesJsonFunc));
	loader.RegisterFunction(ScalarFunction("mol_block_atoms_json", {LogicalType::VARCHAR}, LogicalType::VARCHAR, MolBlockAtomsJsonFunc));
	loader.RegisterFunction(ScalarFunction("mol_block_bonds_json", {LogicalType::VARCHAR}, LogicalType::VARCHAR, MolBlockBondsJsonFunc));
	loader.RegisterFunction(ScalarFunction("mol_block_json",       {LogicalType::VARCHAR}, LogicalType::VARCHAR, MolBlockJsonFunc));
	loader.RegisterFunction(ScalarFunction("mol_block_has_3d",     {LogicalType::VARCHAR}, LogicalType::BOOLEAN, MolBlockHas3dFunc));
	loader.RegisterFunction(ScalarFunction("mol_block_centroid_x", {LogicalType::VARCHAR}, LogicalType::DOUBLE, MolBlockCentroidXFunc));
	loader.RegisterFunction(ScalarFunction("mol_block_centroid_y", {LogicalType::VARCHAR}, LogicalType::DOUBLE, MolBlockCentroidYFunc));
	loader.RegisterFunction(ScalarFunction("mol_block_centroid_z", {LogicalType::VARCHAR}, LogicalType::DOUBLE, MolBlockCentroidZFunc));
	loader.RegisterFunction(ScalarFunction("mol_block_radius_of_gyration", {LogicalType::VARCHAR}, LogicalType::DOUBLE, MolBlockRadiusOfGyrationFunc));
	loader.RegisterFunction(ScalarFunction("mol_block_min_x",      {LogicalType::VARCHAR}, LogicalType::DOUBLE, MolBlockMinXFunc));
	loader.RegisterFunction(ScalarFunction("mol_block_max_x",      {LogicalType::VARCHAR}, LogicalType::DOUBLE, MolBlockMaxXFunc));
	loader.RegisterFunction(ScalarFunction("mol_block_min_y",      {LogicalType::VARCHAR}, LogicalType::DOUBLE, MolBlockMinYFunc));
	loader.RegisterFunction(ScalarFunction("mol_block_max_y",      {LogicalType::VARCHAR}, LogicalType::DOUBLE, MolBlockMaxYFunc));
	loader.RegisterFunction(ScalarFunction("mol_block_min_z",      {LogicalType::VARCHAR}, LogicalType::DOUBLE, MolBlockMinZFunc));
	loader.RegisterFunction(ScalarFunction("mol_block_max_z",      {LogicalType::VARCHAR}, LogicalType::DOUBLE, MolBlockMaxZFunc));
	loader.RegisterFunction(ScalarFunction("sdf_count",            {LogicalType::VARCHAR}, LogicalType::INTEGER, SdfCountFunc));
	loader.RegisterFunction(ScalarFunction("sdf_property",         {LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::VARCHAR}, LogicalType::VARCHAR, SdfPropertyFunc));
	loader.RegisterFunction(ScalarFunction("sdf_properties_json",  {LogicalType::VARCHAR}, LogicalType::VARCHAR, SdfPropertiesJsonFunc));

	// --- PDB/CIF/XYZ structure ---
	loader.RegisterFunction(ScalarFunction("structure_atom_count",    {LogicalType::VARCHAR}, LogicalType::INTEGER, StructureAtomCountFunc));
	loader.RegisterFunction(ScalarFunction("structure_chain_count",   {LogicalType::VARCHAR}, LogicalType::INTEGER, StructureChainCountFunc));
	loader.RegisterFunction(ScalarFunction("structure_residue_count", {LogicalType::VARCHAR}, LogicalType::INTEGER, StructureResidueCountFunc));
	loader.RegisterFunction(ScalarFunction("structure_model_count",   {LogicalType::VARCHAR}, LogicalType::INTEGER, StructureModelCountFunc));
	loader.RegisterFunction(ScalarFunction("structure_centroid_x",    {LogicalType::VARCHAR}, LogicalType::DOUBLE, StructureCentroidXFunc));
	loader.RegisterFunction(ScalarFunction("structure_centroid_y",    {LogicalType::VARCHAR}, LogicalType::DOUBLE, StructureCentroidYFunc));
	loader.RegisterFunction(ScalarFunction("structure_centroid_z",    {LogicalType::VARCHAR}, LogicalType::DOUBLE, StructureCentroidZFunc));
	loader.RegisterFunction(ScalarFunction("structure_radius_of_gyration", {LogicalType::VARCHAR}, LogicalType::DOUBLE, StructureRadiusOfGyrationFunc));
	loader.RegisterFunction(ScalarFunction("structure_min_x",         {LogicalType::VARCHAR}, LogicalType::DOUBLE, StructureMinXFunc));
	loader.RegisterFunction(ScalarFunction("structure_max_x",         {LogicalType::VARCHAR}, LogicalType::DOUBLE, StructureMaxXFunc));
	loader.RegisterFunction(ScalarFunction("structure_min_y",         {LogicalType::VARCHAR}, LogicalType::DOUBLE, StructureMinYFunc));
	loader.RegisterFunction(ScalarFunction("structure_max_y",         {LogicalType::VARCHAR}, LogicalType::DOUBLE, StructureMaxYFunc));
	loader.RegisterFunction(ScalarFunction("structure_min_z",         {LogicalType::VARCHAR}, LogicalType::DOUBLE, StructureMinZFunc));
	loader.RegisterFunction(ScalarFunction("structure_max_z",         {LogicalType::VARCHAR}, LogicalType::DOUBLE, StructureMaxZFunc));

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
