#define DUCKDB_EXTENSION_MAIN

#include "smiles_extension.hpp"
#include "smiles_parser.hpp"

#include "duckdb.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "duckdb/common/vector_operations/unary_executor.hpp"

namespace duckdb {

// ============================================================================
// mol_is_valid(VARCHAR) → BOOLEAN
// ============================================================================
static void MolIsValidFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	UnaryExecutor::Execute<string_t, bool>(args.data[0], result, args.size(), [](string_t input) -> bool {
		std::string smi(input.GetData(), input.GetSize());
		auto mol = smiles::parse(smi);
		return mol.valid;
	});
}

// ============================================================================
// mol_formula(VARCHAR) → VARCHAR  (NULL on invalid)
// ============================================================================
static void MolFormulaFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	UnaryExecutor::ExecuteWithNulls<string_t, string_t>(args.data[0], result, args.size(),
	    [&](string_t input, ValidityMask &mask, idx_t idx) -> string_t {
		    std::string smi(input.GetData(), input.GetSize());
		    auto mol = smiles::parse(smi);
		    if (!mol.valid) {
			    mask.SetInvalid(idx);
			    return string_t();
		    }
		    return StringVector::AddString(result, mol.formula());
	    });
}

// ============================================================================
// mol_num_atoms(VARCHAR) → INTEGER  (heavy atom count, NULL on invalid)
// ============================================================================
static void MolNumAtomsFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	UnaryExecutor::ExecuteWithNulls<string_t, int32_t>(args.data[0], result, args.size(),
	    [](string_t input, ValidityMask &mask, idx_t idx) -> int32_t {
		    std::string smi(input.GetData(), input.GetSize());
		    auto mol = smiles::parse(smi);
		    if (!mol.valid) {
			    mask.SetInvalid(idx);
			    return 0;
		    }
		    return static_cast<int32_t>(mol.heavy_atom_count());
	    });
}

// ============================================================================
// mol_num_bonds(VARCHAR) → INTEGER  (NULL on invalid)
// ============================================================================
static void MolNumBondsFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	UnaryExecutor::ExecuteWithNulls<string_t, int32_t>(args.data[0], result, args.size(),
	    [](string_t input, ValidityMask &mask, idx_t idx) -> int32_t {
		    std::string smi(input.GetData(), input.GetSize());
		    auto mol = smiles::parse(smi);
		    if (!mol.valid) {
			    mask.SetInvalid(idx);
			    return 0;
		    }
		    return static_cast<int32_t>(mol.num_bonds());
	    });
}

// ============================================================================
// Registration
// ============================================================================
void RegisterSmilesFunctions(ExtensionLoader &loader) {
	loader.RegisterFunction(
	    ScalarFunction("mol_is_valid", {LogicalType::VARCHAR}, LogicalType::BOOLEAN, MolIsValidFunction));

	loader.RegisterFunction(
	    ScalarFunction("mol_formula", {LogicalType::VARCHAR}, LogicalType::VARCHAR, MolFormulaFunction));

	loader.RegisterFunction(
	    ScalarFunction("mol_num_atoms", {LogicalType::VARCHAR}, LogicalType::INTEGER, MolNumAtomsFunction));

	loader.RegisterFunction(
	    ScalarFunction("mol_num_bonds", {LogicalType::VARCHAR}, LogicalType::INTEGER, MolNumBondsFunction));
}

// ============================================================================
// Extension lifecycle
// ============================================================================
void SmilesExtension::Load(ExtensionLoader &loader) {
	RegisterSmilesFunctions(loader);
}

std::string SmilesExtension::Name() {
	return "smiles";
}

std::string SmilesExtension::Version() const {
#ifdef EXT_VERSION_SMILES
	return EXT_VERSION_SMILES;
#else
	return "0.1.0";
#endif
}

} // namespace duckdb

extern "C" {

DUCKDB_EXTENSION_API void smiles_duckdb_cpp_init(duckdb::ExtensionLoader &loader) {
	duckdb::SmilesExtension ext;
	ext.Load(loader);
}

DUCKDB_EXTENSION_API void smiles_init(duckdb::DatabaseInstance &db) {
	// Legacy init — redirect handled by DuckDB
}

DUCKDB_EXTENSION_API const char *smiles_version() {
	return duckdb::DuckDB::LibraryVersion();
}
}
