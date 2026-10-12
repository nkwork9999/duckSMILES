#define DUCKDB_EXTENSION_MAIN

#include "ducksmiles_extension.hpp"
#include "ducksmiles.h"
#include "ducksmiles_compat.hpp"

#include "duckdb.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/vector_operations/binary_executor.hpp"
#include "duckdb/common/vector_operations/unary_executor.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "duckdb/parser/parsed_data/create_scalar_function_info.hpp"
#if __has_include("duckdb/common/vector/list_vector.hpp")
#include "duckdb/common/vector/list_vector.hpp"
#else
#include "duckdb/common/types/vector.hpp"
#endif

#include <cmath>
#include <vector>

namespace duckdb {

// ============================================================================
// Helper: call Rust string→bool/int FFI via UnaryExecutor
// ============================================================================

// VARCHAR → BOOLEAN (mol_is_valid, inchi_is_valid, inchikey_is_valid)
#define DEFINE_BOOL_FUNC(FuncName, RustFunc)                                   \
  static void FuncName(DataChunk &args, ExpressionState &state,                \
                       Vector &result) {                                       \
    UnaryExecutor::Execute<string_t, bool>(                                    \
        args.data[0], result, args.size(), [](string_t input) -> bool {        \
          return RustFunc((const uint8_t *)input.GetData(),                    \
                          input.GetSize()) == 1;                               \
        });                                                                    \
  }

// VARCHAR → INTEGER with NULL on -1 sentinel
#define DEFINE_INT_FUNC(FuncName, RustFunc)                                    \
  static void FuncName(DataChunk &args, ExpressionState &state,                \
                       Vector &result) {                                       \
    ducksmiles_compat::ExecuteWithNulls<string_t, int32_t>(                    \
        args.data[0], result, args.size(),                                     \
        [](string_t input, ValidityMask &mask, idx_t idx) -> int32_t {         \
          int32_t val =                                                        \
              RustFunc((const uint8_t *)input.GetData(), input.GetSize());     \
          if (val < 0) {                                                       \
            mask.SetInvalid(idx);                                              \
            return 0;                                                          \
          }                                                                    \
          return val;                                                          \
        });                                                                    \
  }

// VARCHAR → signed INTEGER with NULL on INT32_MIN. Unlike the ordinary count
// sentinel, this preserves legitimate negative formal charges.
#define DEFINE_SIGNED_INT_FUNC(FuncName, RustFunc)                             \
  static void FuncName(DataChunk &args, ExpressionState &state,                \
                       Vector &result) {                                       \
    ducksmiles_compat::ExecuteWithNulls<string_t, int32_t>(                    \
        args.data[0], result, args.size(),                                     \
        [](string_t input, ValidityMask &mask, idx_t idx) -> int32_t {         \
          int32_t val =                                                        \
              RustFunc((const uint8_t *)input.GetData(), input.GetSize());     \
          if (val == INT32_MIN) {                                              \
            mask.SetInvalid(idx);                                              \
            return 0;                                                          \
          }                                                                    \
          return val;                                                          \
        });                                                                    \
  }

// VARCHAR → DOUBLE with NULL on NaN
#define DEFINE_DOUBLE_FUNC(FuncName, RustFunc)                                 \
  static void FuncName(DataChunk &args, ExpressionState &state,                \
                       Vector &result) {                                       \
    ducksmiles_compat::ExecuteWithNulls<string_t, double>(                     \
        args.data[0], result, args.size(),                                     \
        [](string_t input, ValidityMask &mask, idx_t idx) -> double {          \
          double val =                                                         \
              RustFunc((const uint8_t *)input.GetData(), input.GetSize());     \
          if (std::isnan(val)) {                                               \
            mask.SetInvalid(idx);                                              \
            return 0.0;                                                        \
          }                                                                    \
          return val;                                                          \
        });                                                                    \
  }

// VARCHAR → BOOLEAN with NULL on -1 sentinel
#define DEFINE_SENTINEL_BOOL_FUNC(FuncName, RustFunc)                          \
  static void FuncName(DataChunk &args, ExpressionState &state,                \
                       Vector &result) {                                       \
    ducksmiles_compat::ExecuteWithNulls<string_t, bool>(                       \
        args.data[0], result, args.size(),                                     \
        [](string_t input, ValidityMask &mask, idx_t idx) -> bool {            \
          int32_t val =                                                        \
              RustFunc((const uint8_t *)input.GetData(), input.GetSize());     \
          if (val < 0) {                                                       \
            mask.SetInvalid(idx);                                              \
            return false;                                                      \
          }                                                                    \
          return val == 1;                                                     \
        });                                                                    \
  }

// VARCHAR → VARCHAR via Rust buffer-writing FFI. NULL on -1.
#define DEFINE_STR_FUNC(FuncName, RustFunc)                                    \
  static void FuncName(DataChunk &args, ExpressionState &state,                \
                       Vector &result) {                                       \
    ducksmiles_compat::ExecuteWithNulls<string_t, string_t>(                   \
        args.data[0], result, args.size(),                                     \
        [&](string_t input, ValidityMask &mask, idx_t idx) -> string_t {       \
          uint8_t buf[1024];                                                   \
          int32_t len = RustFunc((const uint8_t *)input.GetData(),             \
                                 input.GetSize(), buf, sizeof(buf));           \
          if (len < 0) {                                                       \
            mask.SetInvalid(idx);                                              \
            return string_t();                                                 \
          }                                                                    \
          if (len == 0) {                                                      \
            return StringVector::AddString(result, "");                        \
          }                                                                    \
          if ((size_t)len > sizeof(buf)) {                                     \
            /* The Rust side reports the length it needs when the buffer is    \
               too small and writes nothing, so retry on the heap. */          \
            vector<uint8_t> heap((size_t)len);                                 \
            int32_t written =                                                  \
                RustFunc((const uint8_t *)input.GetData(), input.GetSize(),    \
                         heap.data(), heap.size());                            \
            if (written < 0 || (size_t)written > heap.size()) {                \
              mask.SetInvalid(idx);                                            \
              return string_t();                                               \
            }                                                                  \
            return StringVector::AddString(result, (const char *)heap.data(),  \
                                           written);                           \
          }                                                                    \
          return StringVector::AddString(result, (const char *)buf, len);      \
        });                                                                    \
  }

template <class Callback>
static string_t DynamicStringResult(Vector &result, ValidityMask &mask,
                                    idx_t idx, Callback callback) {
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
    throw InvalidInputException(
        "ducksmiles string FFI length changed between sizing and writing");
  }
  return StringVector::AddString(result, (const char *)buf.data(), buf.size());
}

#define DEFINE_DYNAMIC_STR_FUNC(FuncName, RustFunc)                            \
  static void FuncName(DataChunk &args, ExpressionState &state,                \
                       Vector &result) {                                       \
    ducksmiles_compat::ExecuteWithNulls<string_t, string_t>(                   \
        args.data[0], result, args.size(),                                     \
        [&](string_t input, ValidityMask &mask, idx_t idx) -> string_t {       \
          return DynamicStringResult(                                          \
              result, mask, idx, [&](uint8_t *out, size_t cap) -> int32_t {    \
                return RustFunc((const uint8_t *)input.GetData(),              \
                                input.GetSize(), out, cap);                    \
              });                                                              \
        });                                                                    \
  }

#define DEFINE_BINARY_DYNAMIC_STR_FUNC(FuncName, RustFunc)                     \
  static void FuncName(DataChunk &args, ExpressionState &state,                \
                       Vector &result) {                                       \
    idx_t count = args.size();                                                 \
    ducksmiles_compat::Flatten(args.data[0], count);                           \
    ducksmiles_compat::Flatten(args.data[1], count);                           \
    auto left_data = FlatVector::GetData<string_t>(args.data[0]);              \
    auto right_data = FlatVector::GetData<string_t>(args.data[1]);             \
    auto result_data =                                                         \
        ducksmiles_compat::ResultData<string_t>(result, args.size());          \
    auto &validity = ducksmiles_compat::ResultValidity(result);                \
    auto &left_validity = FlatVector::Validity(args.data[0]);                  \
    auto &right_validity = FlatVector::Validity(args.data[1]);                 \
    for (idx_t i = 0; i < count; i++) {                                        \
      if (!left_validity.RowIsValid(i) || !right_validity.RowIsValid(i)) {     \
        validity.SetInvalid(i);                                                \
        continue;                                                              \
      }                                                                        \
      auto left = left_data[i];                                                \
      auto right = right_data[i];                                              \
      result_data[i] = DynamicStringResult(                                    \
          result, validity, i, [&](uint8_t *out, size_t cap) -> int32_t {      \
            return RustFunc((const uint8_t *)left.GetData(), left.GetSize(),   \
                            (const uint8_t *)right.GetData(), right.GetSize(), \
                            out, cap);                                         \
          });                                                                  \
    }                                                                          \
  }

// ============================================================================
// SMILES functions
// ============================================================================

DEFINE_BOOL_FUNC(MolIsValidFunc, ds_mol_is_valid)
DEFINE_INT_FUNC(MolNumAtomsFunc, ds_mol_num_atoms)
DEFINE_INT_FUNC(MolNumFragmentsFunc, ds_mol_num_fragments)
DEFINE_INT_FUNC(MolNumBondsFunc, ds_mol_num_bonds)
DEFINE_SIGNED_INT_FUNC(MolFormalChargeFunc, ds_mol_formal_charge)
DEFINE_INT_FUNC(MolNumExplicitHFunc, ds_mol_num_explicit_h)
DEFINE_INT_FUNC(MolNumImplicitHFunc, ds_mol_num_implicit_h)
DEFINE_INT_FUNC(MolNumTotalHFunc, ds_mol_num_total_h)
DEFINE_INT_FUNC(MolNumSingleBondsFunc, ds_mol_num_single_bonds)
DEFINE_INT_FUNC(MolNumDoubleBondsFunc, ds_mol_num_double_bonds)
DEFINE_INT_FUNC(MolNumTripleBondsFunc, ds_mol_num_triple_bonds)
DEFINE_INT_FUNC(MolNumAromaticBondsFunc, ds_mol_num_aromatic_bonds)
DEFINE_INT_FUNC(MolNumRingAtomsFunc, ds_mol_num_ring_atoms)
DEFINE_INT_FUNC(MolNumRingBondsFunc, ds_mol_num_ring_bonds)
DEFINE_INT_FUNC(MolLargestRingSizeFunc, ds_mol_largest_ring_size)
DEFINE_INT_FUNC(MolNumAromaticAtomsFunc, ds_mol_num_aromatic_atoms)
DEFINE_INT_FUNC(MolNumCarbonsFunc, ds_mol_num_carbons)
DEFINE_INT_FUNC(MolNumNitrogensFunc, ds_mol_num_nitrogens)
DEFINE_INT_FUNC(MolNumOxygensFunc, ds_mol_num_oxygens)
DEFINE_INT_FUNC(MolNumHalogensFunc, ds_mol_num_halogens)
DEFINE_DOUBLE_FUNC(MolHeteroatomFractionFunc, ds_mol_heteroatom_fraction)
DEFINE_DOUBLE_FUNC(MolAromaticFractionFunc, ds_mol_aromatic_fraction)
DEFINE_DOUBLE_FUNC(MolHeavyAtomMassFunc, ds_mol_heavy_atom_mass)
DEFINE_DOUBLE_FUNC(MolMeanDegreeFunc, ds_mol_mean_degree)
DEFINE_STR_FUNC(MolFormulaFunc, ds_mol_formula)
DEFINE_DOUBLE_FUNC(MolWeightFunc, ds_mol_weight)
DEFINE_DOUBLE_FUNC(MolExactMassFunc, ds_mol_exact_mass)
DEFINE_DOUBLE_FUNC(LogpCrippenFunc, ds_logp_crippen)
DEFINE_DOUBLE_FUNC(TpsaFunc, ds_tpsa)
DEFINE_STR_FUNC(CanonicalSmilesFunc, ds_canonical_smiles)
DEFINE_DYNAMIC_STR_FUNC(MurckoScaffoldFunc, ds_murcko_scaffold)
DEFINE_DYNAMIC_STR_FUNC(GenericScaffoldFunc, ds_generic_scaffold)
DEFINE_DYNAMIC_STR_FUNC(RingSystemsJsonFunc, ds_ring_systems_json)
DEFINE_BINARY_DYNAMIC_STR_FUNC(MolHashFunc, ds_mol_hash)
DEFINE_DYNAMIC_STR_FUNC(LargestFragmentFunc, ds_largest_fragment)
DEFINE_DYNAMIC_STR_FUNC(StripSaltsFunc, ds_strip_salts)
DEFINE_DYNAMIC_STR_FUNC(NeutralizeChargesFunc, ds_neutralize_charges)
DEFINE_DYNAMIC_STR_FUNC(NormalizeSmilesFunc, ds_normalize_smiles)
DEFINE_DYNAMIC_STR_FUNC(FragmentParentFunc, ds_fragment_parent)
DEFINE_BINARY_DYNAMIC_STR_FUNC(McsSmartsFunc, ds_mcs_smarts)
DEFINE_BINARY_DYNAMIC_STR_FUNC(McsJsonFunc, ds_mcs_json)
DEFINE_DYNAMIC_STR_FUNC(ScaffoldNetworkJsonFunc, ds_scaffold_network_json)
DEFINE_INT_FUNC(NumHAcceptorsFunc, ds_num_h_acceptors)
DEFINE_INT_FUNC(NumHDonorsFunc, ds_num_h_donors)
DEFINE_INT_FUNC(NumRotatableBondsFunc, ds_num_rotatable_bonds)
DEFINE_INT_FUNC(RingCountFunc, ds_ring_count)
DEFINE_INT_FUNC(NumAromaticRingsFunc, ds_num_aromatic_rings)
DEFINE_INT_FUNC(NumHeteroatomsFunc, ds_num_heteroatoms)
DEFINE_DOUBLE_FUNC(FractionCsp3Func, ds_fraction_csp3)
DEFINE_DOUBLE_FUNC(MolMrFunc, ds_mol_mr)
DEFINE_DOUBLE_FUNC(QedFunc, ds_qed)
DEFINE_INT_FUNC(NumAliphaticRingsFunc, ds_num_aliphatic_rings)
DEFINE_INT_FUNC(NumSaturatedRingsFunc, ds_num_saturated_rings)
DEFINE_INT_FUNC(NumAromaticHeterocyclesFunc, ds_num_aromatic_heterocycles)
DEFINE_INT_FUNC(NumAromaticCarbocyclesFunc, ds_num_aromatic_carbocycles)
DEFINE_INT_FUNC(NumSaturatedHeterocyclesFunc, ds_num_saturated_heterocycles)
DEFINE_INT_FUNC(NumSaturatedCarbocyclesFunc, ds_num_saturated_carbocycles)
DEFINE_INT_FUNC(NumAliphaticHeterocyclesFunc, ds_num_aliphatic_heterocycles)
DEFINE_INT_FUNC(NumAliphaticCarbocyclesFunc, ds_num_aliphatic_carbocycles)

// ── ADMET / drug-likeness rule panels + toxicophore structural alerts ───────
DEFINE_DYNAMIC_STR_FUNC(AdmetJsonFunc, ds_admet_json)
DEFINE_DYNAMIC_STR_FUNC(StructuralAlertsJsonFunc, ds_structural_alerts_json)
DEFINE_INT_FUNC(StructuralAlertCountFunc, ds_structural_alert_count)
DEFINE_INT_FUNC(LipinskiViolationsFunc, ds_lipinski_violations)
// Protein PDB → PDBQT (Vina atom typing).
DEFINE_DYNAMIC_STR_FUNC(PdbToPdbqtFunc, ds_pdb_to_pdbqt)

// druglikeness_pass(smiles, rule) → INTEGER (1 pass / 0 fail / -1 invalid /
// -2 unknown rule). Binary VARCHAR inputs, integer output.
static void DruglikenessPassFunc(DataChunk &args, ExpressionState &state,
                                 Vector &result) {
  BinaryExecutor::Execute<string_t, string_t, int32_t>(
      args.data[0], args.data[1], result, args.size(),
      [](string_t smi, string_t rule) -> int32_t {
        return ds_druglikeness_pass(
            (const uint8_t *)smi.GetData(), smi.GetSize(),
            (const uint8_t *)rule.GetData(), rule.GetSize());
      });
}

// ── Docking pipeline (conformer / PDBQT / dock) ─────────────────────────────

// smiles_to_pdbqt(smiles, seed) → VARCHAR (ligand PDBQT). NULL on invalid.
static void SmilesToPdbqtFunc(DataChunk &args, ExpressionState &state,
                              Vector &result) {
  idx_t count = args.size();
  ducksmiles_compat::Flatten(args.data[0], count);
  ducksmiles_compat::Flatten(args.data[1], count);
  auto smi = FlatVector::GetData<string_t>(args.data[0]);
  auto seed = FlatVector::GetData<int64_t>(args.data[1]);
  auto &smi_valid = FlatVector::Validity(args.data[0]);
  auto out = ducksmiles_compat::ResultData<string_t>(result, args.size());
  auto &out_valid = ducksmiles_compat::ResultValidity(result);
  std::vector<uint8_t> buf(1u << 20); // 1 MiB
  for (idx_t i = 0; i < count; i++) {
    if (!smi_valid.RowIsValid(i)) {
      out_valid.SetInvalid(i);
      continue;
    }
    int32_t n =
        ds_smiles_to_pdbqt((const uint8_t *)smi[i].GetData(), smi[i].GetSize(),
                           (uint64_t)seed[i], buf.data(), buf.size());
    if (n < 0) {
      out_valid.SetInvalid(i);
      continue;
    }
    out[i] =
        StringVector::AddString(result, (const char *)buf.data(), (size_t)n);
  }
}

// dock(smiles, pdb, cx,cy,cz, sx,sy,sz, n_runs, seed [, ph]) → VARCHAR JSON.
// Accepts 10 args (ph defaults to 7.4) or 11 args (explicit ph).
static void DockFunc(DataChunk &args, ExpressionState &state, Vector &result) {
  idx_t count = args.size();
  idx_t ncol = args.ColumnCount();
  for (idx_t c = 0; c < ncol; c++)
    ducksmiles_compat::Flatten(args.data[c], count);
  auto smi = FlatVector::GetData<string_t>(args.data[0]);
  auto pdb = FlatVector::GetData<string_t>(args.data[1]);
  auto cx = FlatVector::GetData<double>(args.data[2]);
  auto cy = FlatVector::GetData<double>(args.data[3]);
  auto cz = FlatVector::GetData<double>(args.data[4]);
  auto sx = FlatVector::GetData<double>(args.data[5]);
  auto sy = FlatVector::GetData<double>(args.data[6]);
  auto sz = FlatVector::GetData<double>(args.data[7]);
  auto nruns = FlatVector::GetData<int32_t>(args.data[8]);
  auto seed = FlatVector::GetData<int64_t>(args.data[9]);
  const double *phcol =
      (ncol >= 11) ? FlatVector::GetData<double>(args.data[10]) : nullptr;
  auto &smi_valid = FlatVector::Validity(args.data[0]);
  auto &pdb_valid = FlatVector::Validity(args.data[1]);
  auto out = ducksmiles_compat::ResultData<string_t>(result, args.size());
  auto &out_valid = ducksmiles_compat::ResultValidity(result);
  std::vector<uint8_t> buf(1u << 20); // 1 MiB
  for (idx_t i = 0; i < count; i++) {
    if (!smi_valid.RowIsValid(i) || !pdb_valid.RowIsValid(i)) {
      out_valid.SetInvalid(i);
      continue;
    }
    double ph = phcol ? phcol[i] : 7.4;
    int32_t n =
        ds_dock((const uint8_t *)smi[i].GetData(), smi[i].GetSize(),
                (const uint8_t *)pdb[i].GetData(), pdb[i].GetSize(), cx[i],
                cy[i], cz[i], sx[i], sy[i], sz[i], (uint32_t)nruns[i],
                (uint64_t)seed[i], ph, buf.data(), buf.size());
    if (n < 0) {
      out_valid.SetInvalid(i);
      continue;
    }
    out[i] =
        StringVector::AddString(result, (const char *)buf.data(), (size_t)n);
  }
}

// ── Virtual-screening benchmark metrics over LIST<DOUBLE>, LIST<BOOLEAN> ────
// roc_auc(scores, labels), enrichment_factor(scores, labels, fraction),
// bedroc(scores, labels, alpha). One DOUBLE per row; NaN on invalid input.
// Reads each row's two lists into temp arrays and calls the Rust core.

template <class Metric>
static void BenchmarkMetric(DataChunk &args, Vector &result, Metric metric) {
  idx_t count = args.size();
  auto &scores_list = args.data[0];
  auto &labels_list = args.data[1];

  UnifiedVectorFormat s_fmt, l_fmt;
  ducksmiles_compat::ToUnifiedFormat(scores_list, count, s_fmt);
  ducksmiles_compat::ToUnifiedFormat(labels_list, count, l_fmt);
  auto s_entries = UnifiedVectorFormat::GetData<list_entry_t>(s_fmt);
  auto l_entries = UnifiedVectorFormat::GetData<list_entry_t>(l_fmt);

  // Flatten child vectors so we can index them directly.
  auto &s_child = ducksmiles_compat::ListChild(scores_list);
  auto &l_child = ducksmiles_compat::ListChild(labels_list);
  idx_t s_child_len = ListVector::GetListSize(scores_list);
  idx_t l_child_len = ListVector::GetListSize(labels_list);
  ducksmiles_compat::Flatten(s_child, s_child_len);
  ducksmiles_compat::Flatten(l_child, l_child_len);
  auto s_data = FlatVector::GetData<double>(s_child);
  auto l_data = FlatVector::GetData<bool>(l_child);

  auto out = ducksmiles_compat::ResultData<double>(result, args.size());
  auto &out_valid = ducksmiles_compat::ResultValidity(result);

  for (idx_t i = 0; i < count; i++) {
    auto si = s_fmt.sel->get_index(i);
    auto li = l_fmt.sel->get_index(i);
    if (!s_fmt.validity.RowIsValid(si) || !l_fmt.validity.RowIsValid(li)) {
      out_valid.SetInvalid(i);
      continue;
    }
    auto se = s_entries[si];
    auto le = l_entries[li];
    idx_t n = se.length;
    if (n == 0 || le.length != n) {
      out_valid.SetInvalid(i);
      continue;
    }
    std::vector<double> scores(n);
    std::vector<uint8_t> labels(n);
    for (idx_t k = 0; k < n; k++) {
      scores[k] = s_data[se.offset + k];
      labels[k] = l_data[le.offset + k] ? 1 : 0;
    }
    out[i] = metric(scores.data(), labels.data(), n, i);
  }
}

static void RocAucFunc(DataChunk &args, ExpressionState &state,
                       Vector &result) {
  BenchmarkMetric(args, result,
                  [](const double *s, const uint8_t *l, size_t n, idx_t) {
                    return ds_roc_auc(s, l, n);
                  });
}

static void EnrichmentFactorFunc(DataChunk &args, ExpressionState &state,
                                 Vector &result) {
  idx_t count = args.size();
  ducksmiles_compat::Flatten(args.data[2], count);
  auto frac = FlatVector::GetData<double>(args.data[2]);
  BenchmarkMetric(args, result,
                  [&](const double *s, const uint8_t *l, size_t n, idx_t i) {
                    return ds_enrichment_factor(s, l, n, frac[i]);
                  });
}

static void BedrocFunc(DataChunk &args, ExpressionState &state,
                       Vector &result) {
  idx_t count = args.size();
  ducksmiles_compat::Flatten(args.data[2], count);
  auto alpha = FlatVector::GetData<double>(args.data[2]);
  BenchmarkMetric(args, result,
                  [&](const double *s, const uint8_t *l, size_t n, idx_t i) {
                    return ds_bedroc(s, l, n, alpha[i]);
                  });
}

// prepare_receptor(pdb, ph) → VARCHAR PDBQT (protonated + polar-H added).
static void PrepareReceptorFunc(DataChunk &args, ExpressionState &state,
                                Vector &result) {
  idx_t count = args.size();
  ducksmiles_compat::Flatten(args.data[0], count);
  ducksmiles_compat::Flatten(args.data[1], count);
  auto pdb = FlatVector::GetData<string_t>(args.data[0]);
  auto ph = FlatVector::GetData<double>(args.data[1]);
  auto &pdb_valid = FlatVector::Validity(args.data[0]);
  auto out = ducksmiles_compat::ResultData<string_t>(result, args.size());
  auto &out_valid = ducksmiles_compat::ResultValidity(result);
  for (idx_t i = 0; i < count; i++) {
    if (!pdb_valid.RowIsValid(i)) {
      out_valid.SetInvalid(i);
      continue;
    }
    out[i] = DynamicStringResult(
        result, out_valid, i, [&](uint8_t *o, size_t cap) -> int32_t {
          return ds_prepare_receptor((const uint8_t *)pdb[i].GetData(),
                                     pdb[i].GetSize(), ph[i], o, cap);
        });
  }
}

static void MolHashMethodsJsonFunc(DataChunk &args, ExpressionState &state,
                                   Vector &result) {
  idx_t count = args.size();
  auto result_data =
      ducksmiles_compat::ResultData<string_t>(result, args.size());
  auto &validity = ducksmiles_compat::ResultValidity(result);
  for (idx_t i = 0; i < count; i++) {
    result_data[i] = DynamicStringResult(
        result, validity, i, [&](uint8_t *out, size_t cap) -> int32_t {
          return ds_mol_hash_methods_json(out, cap);
        });
  }
}

static void MolHasSubstructureFunc(DataChunk &args, ExpressionState &state,
                                   Vector &result) {
  idx_t count = args.size();
  ducksmiles_compat::Flatten(args.data[0], count);
  ducksmiles_compat::Flatten(args.data[1], count);
  auto smi_data = FlatVector::GetData<string_t>(args.data[0]);
  auto smarts_data = FlatVector::GetData<string_t>(args.data[1]);
  auto result_data = ducksmiles_compat::ResultData<bool>(result, args.size());
  auto &validity = ducksmiles_compat::ResultValidity(result);
  for (idx_t i = 0; i < count; i++) {
    if (!FlatVector::Validity(args.data[0]).RowIsValid(i) ||
        !FlatVector::Validity(args.data[1]).RowIsValid(i)) {
      validity.SetInvalid(i);
      continue;
    }
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

static void MolSubstructureCountFunc(DataChunk &args, ExpressionState &state,
                                     Vector &result) {
  idx_t count = args.size();
  ducksmiles_compat::Flatten(args.data[0], count);
  ducksmiles_compat::Flatten(args.data[1], count);
  auto smi_data = FlatVector::GetData<string_t>(args.data[0]);
  auto smarts_data = FlatVector::GetData<string_t>(args.data[1]);
  auto result_data =
      ducksmiles_compat::ResultData<int32_t>(result, args.size());
  auto &validity = ducksmiles_compat::ResultValidity(result);
  for (idx_t i = 0; i < count; i++) {
    if (!FlatVector::Validity(args.data[0]).RowIsValid(i) ||
        !FlatVector::Validity(args.data[1]).RowIsValid(i)) {
      validity.SetInvalid(i);
      continue;
    }
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

static void MolSubstructureMatchesJsonFunc(DataChunk &args,
                                           ExpressionState &state,
                                           Vector &result) {
  idx_t count = args.size();
  ducksmiles_compat::Flatten(args.data[0], count);
  ducksmiles_compat::Flatten(args.data[1], count);
  auto smi_data = FlatVector::GetData<string_t>(args.data[0]);
  auto smarts_data = FlatVector::GetData<string_t>(args.data[1]);
  auto result_data =
      ducksmiles_compat::ResultData<string_t>(result, args.size());
  auto &validity = ducksmiles_compat::ResultValidity(result);
  auto &smi_validity = FlatVector::Validity(args.data[0]);
  auto &smarts_validity = FlatVector::Validity(args.data[1]);
  for (idx_t i = 0; i < count; i++) {
    if (!smi_validity.RowIsValid(i) || !smarts_validity.RowIsValid(i)) {
      validity.SetInvalid(i);
      continue;
    }
    auto smi = smi_data[i];
    auto smarts = smarts_data[i];
    result_data[i] = DynamicStringResult(
        result, validity, i, [&](uint8_t *out, size_t cap) -> int32_t {
          return ds_mol_substructure_matches_json(
              (const uint8_t *)smi.GetData(), smi.GetSize(),
              (const uint8_t *)smarts.GetData(), smarts.GetSize(), out, cap);
        });
  }
}

// add_hydrogens uses a larger 16KB buffer to handle drug-sized molecules
// (verbose SMILES with all H broken out can be ~5x the heavy-atom SMILES
// length).
static void AddHydrogensFunc(DataChunk &args, ExpressionState &state,
                             Vector &result) {
  ducksmiles_compat::ExecuteWithNulls<string_t, string_t>(
      args.data[0], result, args.size(),
      [&](string_t input, ValidityMask &mask, idx_t idx) -> string_t {
        uint8_t buf[16384];
        int32_t len = ds_add_hydrogens((const uint8_t *)input.GetData(),
                                       input.GetSize(), buf, sizeof(buf));
        if (len < 0) {
          mask.SetInvalid(idx);
          return string_t();
        }
        if (len == 0) {
          return StringVector::AddString(result, "");
        }
        if ((size_t)len > sizeof(buf)) {
          // Buffer too small: the Rust side reported the size it needs.
          vector<uint8_t> heap((size_t)len);
          int32_t written =
              ds_add_hydrogens((const uint8_t *)input.GetData(),
                               input.GetSize(), heap.data(), heap.size());
          if (written < 0 || (size_t)written > heap.size()) {
            mask.SetInvalid(idx);
            return string_t();
          }
          return StringVector::AddString(result, (const char *)heap.data(),
                                         written);
        }
        return StringVector::AddString(result, (const char *)buf, len);
      });
}

// morgan_fp_bits(smi, radius, n_bits) → BLOB (ceil(n_bits/8) bytes).
// 16 KiB buffer is enough for up to 131072-bit fingerprints (standard sizes are
// 1024–4096).
static constexpr size_t MORGAN_BUF_BYTES = 16384;

static void MorganFpBitsFunc3(DataChunk &args, ExpressionState &state,
                              Vector &result) {
  idx_t count = args.size();
  ducksmiles_compat::Flatten(args.data[0], count);
  ducksmiles_compat::Flatten(args.data[1], count);
  ducksmiles_compat::Flatten(args.data[2], count);
  auto smi_data = FlatVector::GetData<string_t>(args.data[0]);
  auto radius_data = FlatVector::GetData<int32_t>(args.data[1]);
  auto nbits_data = FlatVector::GetData<int32_t>(args.data[2]);
  auto out = ducksmiles_compat::ResultData<string_t>(result, count);
  auto &validity = ducksmiles_compat::ResultValidity(result);
  for (idx_t i = 0; i < count; i++) {
    if (!FlatVector::Validity(args.data[0]).RowIsValid(i) ||
        !FlatVector::Validity(args.data[1]).RowIsValid(i) ||
        !FlatVector::Validity(args.data[2]).RowIsValid(i)) {
      validity.SetInvalid(i);
      continue;
    }
    int32_t r = radius_data[i];
    int32_t n = nbits_data[i];
    if (r < 0 || n <= 0) {
      validity.SetInvalid(i);
      continue;
    }
    uint8_t buf[MORGAN_BUF_BYTES];
    int32_t len = ds_morgan_fp_bits((const uint8_t *)smi_data[i].GetData(),
                                    smi_data[i].GetSize(), (uint32_t)r,
                                    (uint32_t)n, buf, sizeof(buf));
    if (len < 0) {
      validity.SetInvalid(i);
      continue;
    }
    auto blob = StringVector::AddStringOrBlob(result, (const char *)buf, len);
    out[i] = blob;
  }
}

// morgan_fp_bits(smi) → BLOB with defaults: ECFP4 (radius=2, 2048 bits → 256
// bytes).
static void MorganFpBitsFunc1(DataChunk &args, ExpressionState &state,
                              Vector &result) {
  ducksmiles_compat::ExecuteWithNulls<string_t, string_t>(
      args.data[0], result, args.size(),
      [&](string_t input, ValidityMask &mask, idx_t idx) -> string_t {
        uint8_t buf[MORGAN_BUF_BYTES];
        int32_t len =
            ds_morgan_fp_bits((const uint8_t *)input.GetData(), input.GetSize(),
                              2u, 2048u, buf, sizeof(buf));
        if (len < 0) {
          mask.SetInvalid(idx);
          return string_t();
        }
        return StringVector::AddStringOrBlob(result, (const char *)buf, len);
      });
}

// maccs_keys(smi) → BLOB (fixed 21 bytes, 167-bit MACCS vector).
static constexpr size_t MACCS_BUF_BYTES = 21;

static void MaccsKeysFunc(DataChunk &args, ExpressionState &state,
                          Vector &result) {
  ducksmiles_compat::ExecuteWithNulls<string_t, string_t>(
      args.data[0], result, args.size(),
      [&](string_t input, ValidityMask &mask, idx_t idx) -> string_t {
        uint8_t buf[MACCS_BUF_BYTES];
        int32_t len = ds_maccs_keys((const uint8_t *)input.GetData(),
                                    input.GetSize(), buf, sizeof(buf));
        if (len < 0) {
          mask.SetInvalid(idx);
          return string_t();
        }
        return StringVector::AddStringOrBlob(result, (const char *)buf, len);
      });
}

// tanimoto_bit(BLOB, BLOB) → DOUBLE.
// Computes popcount(a & b) / popcount(a | b) on the raw BLOB bytes —
// no CAST AS BIT, no per-row table dispatch. Both arguments must be the
// same length (typically the 256-byte BLOB produced by morgan_fp_bits).
// BinaryExecutor handles flat / constant / dictionary vectors and NULL
// propagation; we only worry about length-mismatch (clear user error).
static void TanimotoBitFunc(DataChunk &args, ExpressionState &state,
                            Vector &result) {
  BinaryExecutor::Execute<string_t, string_t, double>(
      args.data[0], args.data[1], result, args.size(),
      [](string_t a, string_t b) -> double {
        if (a.GetSize() != b.GetSize()) {
          throw InvalidInputException(
              "tanimoto_bit: BLOB lengths differ (%llu vs %llu bytes)",
              (unsigned long long)a.GetSize(), (unsigned long long)b.GetSize());
        }
        return ds_tanimoto_bit((const uint8_t *)a.GetData(), a.GetSize(),
                               (const uint8_t *)b.GetData(), b.GetSize());
      });
}

// <name>_bit(BLOB, BLOB) → DOUBLE for the rest of the symmetric similarity
// family. Same length-mismatch contract as tanimoto_bit; SqlName is the label
// shown in the error so the user sees the function they actually called.
#define DEFINE_SIMILARITY_FUNC(FuncName, RustFunc, SqlName)                    \
  static void FuncName(DataChunk &args, ExpressionState &state,                \
                       Vector &result) {                                       \
    BinaryExecutor::Execute<string_t, string_t, double>(                       \
        args.data[0], args.data[1], result, args.size(),                       \
        [](string_t a, string_t b) -> double {                                 \
          if (a.GetSize() != b.GetSize()) {                                    \
            throw InvalidInputException(                                       \
                SqlName ": BLOB lengths differ (%llu vs %llu bytes)",          \
                (unsigned long long)a.GetSize(),                               \
                (unsigned long long)b.GetSize());                              \
          }                                                                    \
          return RustFunc((const uint8_t *)a.GetData(), a.GetSize(),           \
                          (const uint8_t *)b.GetData(), b.GetSize());          \
        });                                                                    \
  }

DEFINE_SIMILARITY_FUNC(DiceBitFunc, ds_dice_bit, "dice_bit")
DEFINE_SIMILARITY_FUNC(CosineBitFunc, ds_cosine_bit, "cosine_bit")
DEFINE_SIMILARITY_FUNC(KulczynskiBitFunc, ds_kulczynski_bit, "kulczynski_bit")
DEFINE_SIMILARITY_FUNC(SokalBitFunc, ds_sokal_bit, "sokal_bit")
DEFINE_SIMILARITY_FUNC(McConnaugheyBitFunc, ds_mcconnaughey_bit,
                       "mcconnaughey_bit")
DEFINE_SIMILARITY_FUNC(AsymmetricBitFunc, ds_asymmetric_bit, "asymmetric_bit")
DEFINE_SIMILARITY_FUNC(BraunBlanquetBitFunc, ds_braun_blanquet_bit,
                       "braun_blanquet_bit")
DEFINE_SIMILARITY_FUNC(RusselBitFunc, ds_russel_bit, "russel_bit")

// tversky_bit(BLOB, BLOB, alpha, beta) → DOUBLE. alpha/beta weight the two
// fingerprints (both in [0, 1]); alpha=beta=1 reduces to Tanimoto, 0.5/0.5 to
// Dice. DuckDB has no 4-ary executor, so we read each argument via
// UnifiedVectorFormat (handles flat / constant / dictionary vectors + NULLs).
static void TverskyBitFunc(DataChunk &args, ExpressionState &state,
                           Vector &result) {
  auto count = args.size();
  UnifiedVectorFormat a_fmt, b_fmt, alpha_fmt, beta_fmt;
  ducksmiles_compat::ToUnifiedFormat(args.data[0], count, a_fmt);
  ducksmiles_compat::ToUnifiedFormat(args.data[1], count, b_fmt);
  ducksmiles_compat::ToUnifiedFormat(args.data[2], count, alpha_fmt);
  ducksmiles_compat::ToUnifiedFormat(args.data[3], count, beta_fmt);

  auto a_vals = UnifiedVectorFormat::GetData<string_t>(a_fmt);
  auto b_vals = UnifiedVectorFormat::GetData<string_t>(b_fmt);
  auto alpha_vals = UnifiedVectorFormat::GetData<double>(alpha_fmt);
  auto beta_vals = UnifiedVectorFormat::GetData<double>(beta_fmt);

  result.SetVectorType(VectorType::FLAT_VECTOR);
  auto out = ducksmiles_compat::ResultData<double>(result, args.size());
  auto &out_validity = ducksmiles_compat::ResultValidity(result);

  for (idx_t i = 0; i < count; i++) {
    auto ai = a_fmt.sel->get_index(i);
    auto bi = b_fmt.sel->get_index(i);
    auto pi = alpha_fmt.sel->get_index(i);
    auto qi = beta_fmt.sel->get_index(i);
    if (!a_fmt.validity.RowIsValid(ai) || !b_fmt.validity.RowIsValid(bi) ||
        !alpha_fmt.validity.RowIsValid(pi) ||
        !beta_fmt.validity.RowIsValid(qi)) {
      out_validity.SetInvalid(i);
      continue;
    }
    auto &a = a_vals[ai];
    auto &b = b_vals[bi];
    if (a.GetSize() != b.GetSize()) {
      throw InvalidInputException(
          "tversky_bit: BLOB lengths differ (%llu vs %llu bytes)",
          (unsigned long long)a.GetSize(), (unsigned long long)b.GetSize());
    }
    double alpha = alpha_vals[pi];
    double beta = beta_vals[qi];
    double sim =
        ds_tversky_bit((const uint8_t *)a.GetData(), a.GetSize(),
                       (const uint8_t *)b.GetData(), b.GetSize(), alpha, beta);
    if (std::isnan(sim)) {
      throw InvalidInputException(
          "tversky_bit: alpha and beta must both be in [0, 1] "
          "(got alpha=%g, beta=%g)",
          alpha, beta);
    }
    out[i] = sim;
  }
  if (count == 1) {
    result.SetVectorType(VectorType::CONSTANT_VECTOR);
  }
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
static void InchiSkeletonMatchFunc(DataChunk &args, ExpressionState &state,
                                   Vector &result) {
  idx_t count = args.size();
  ducksmiles_compat::Flatten(args.data[0], count);
  ducksmiles_compat::Flatten(args.data[1], count);
  auto a_data = FlatVector::GetData<string_t>(args.data[0]);
  auto b_data = FlatVector::GetData<string_t>(args.data[1]);
  auto result_data = ducksmiles_compat::ResultData<bool>(result, args.size());
  auto &validity = ducksmiles_compat::ResultValidity(result);
  for (idx_t i = 0; i < count; i++) {
    if (!FlatVector::Validity(args.data[0]).RowIsValid(i) ||
        !FlatVector::Validity(args.data[1]).RowIsValid(i)) {
      validity.SetInvalid(i);
      continue;
    }
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
DEFINE_DYNAMIC_STR_FUNC(MolBlockPropertiesJsonFunc,
                        ds_mol_block_properties_json)
DEFINE_DYNAMIC_STR_FUNC(MolBlockAtomsJsonFunc, ds_mol_block_atoms_json)
DEFINE_DYNAMIC_STR_FUNC(MolBlockBondsJsonFunc, ds_mol_block_bonds_json)
DEFINE_DYNAMIC_STR_FUNC(MolBlockJsonFunc, ds_mol_block_json)
DEFINE_DYNAMIC_STR_FUNC(SdfPropertiesJsonFunc, ds_sdf_properties_json)
DEFINE_INT_FUNC(SdfCountFunc, ds_sdf_count)
DEFINE_SENTINEL_BOOL_FUNC(MolBlockHas3dFunc, ds_mol_block_has_3d)
DEFINE_DOUBLE_FUNC(MolBlockCentroidXFunc, ds_mol_block_centroid_x)
DEFINE_DOUBLE_FUNC(MolBlockCentroidYFunc, ds_mol_block_centroid_y)
DEFINE_DOUBLE_FUNC(MolBlockCentroidZFunc, ds_mol_block_centroid_z)
DEFINE_DOUBLE_FUNC(MolBlockRadiusOfGyrationFunc,
                   ds_mol_block_radius_of_gyration)
DEFINE_DOUBLE_FUNC(MolBlockMinXFunc, ds_mol_block_min_x)
DEFINE_DOUBLE_FUNC(MolBlockMaxXFunc, ds_mol_block_max_x)
DEFINE_DOUBLE_FUNC(MolBlockMinYFunc, ds_mol_block_min_y)
DEFINE_DOUBLE_FUNC(MolBlockMaxYFunc, ds_mol_block_max_y)
DEFINE_DOUBLE_FUNC(MolBlockMinZFunc, ds_mol_block_min_z)
DEFINE_DOUBLE_FUNC(MolBlockMaxZFunc, ds_mol_block_max_z)

static void MolBlockPropertyFunc(DataChunk &args, ExpressionState &state,
                                 Vector &result) {
  idx_t count = args.size();
  ducksmiles_compat::Flatten(args.data[0], count);
  ducksmiles_compat::Flatten(args.data[1], count);
  auto mol_data = FlatVector::GetData<string_t>(args.data[0]);
  auto key_data = FlatVector::GetData<string_t>(args.data[1]);
  auto result_data =
      ducksmiles_compat::ResultData<string_t>(result, args.size());
  auto &validity = ducksmiles_compat::ResultValidity(result);
  auto &mol_validity = FlatVector::Validity(args.data[0]);
  auto &key_validity = FlatVector::Validity(args.data[1]);
  for (idx_t i = 0; i < count; i++) {
    if (!mol_validity.RowIsValid(i) || !key_validity.RowIsValid(i)) {
      validity.SetInvalid(i);
      continue;
    }
    auto mol = mol_data[i];
    auto key = key_data[i];
    result_data[i] = DynamicStringResult(
        result, validity, i, [&](uint8_t *out, size_t cap) -> int32_t {
          return ds_mol_block_property(
              (const uint8_t *)mol.GetData(), mol.GetSize(),
              (const uint8_t *)key.GetData(), key.GetSize(), out, cap);
        });
  }
}

static void SdfPropertyFunc(DataChunk &args, ExpressionState &state,
                            Vector &result) {
  idx_t count = args.size();
  ducksmiles_compat::Flatten(args.data[0], count);
  ducksmiles_compat::Flatten(args.data[1], count);
  ducksmiles_compat::Flatten(args.data[2], count);
  auto sdf_data = FlatVector::GetData<string_t>(args.data[0]);
  auto index_data = FlatVector::GetData<int32_t>(args.data[1]);
  auto key_data = FlatVector::GetData<string_t>(args.data[2]);
  auto result_data =
      ducksmiles_compat::ResultData<string_t>(result, args.size());
  auto &validity = ducksmiles_compat::ResultValidity(result);
  auto &sdf_validity = FlatVector::Validity(args.data[0]);
  auto &index_validity = FlatVector::Validity(args.data[1]);
  auto &key_validity = FlatVector::Validity(args.data[2]);
  for (idx_t i = 0; i < count; i++) {
    if (!sdf_validity.RowIsValid(i) || !index_validity.RowIsValid(i) ||
        !key_validity.RowIsValid(i)) {
      validity.SetInvalid(i);
      continue;
    }
    auto sdf = sdf_data[i];
    auto key = key_data[i];
    result_data[i] = DynamicStringResult(
        result, validity, i, [&](uint8_t *out, size_t cap) -> int32_t {
          return ds_sdf_property((const uint8_t *)sdf.GetData(), sdf.GetSize(),
                                 index_data[i], (const uint8_t *)key.GetData(),
                                 key.GetSize(), out, cap);
        });
  }
}

// ============================================================================
// PDB/CIF/XYZ functions — format arg hardcoded to 0 (auto-detect)
// ============================================================================

#define DEFINE_STRUCTURE_INT_FUNC(FuncName, RustFunc)                          \
  static void FuncName(DataChunk &args, ExpressionState &state,                \
                       Vector &result) {                                       \
    ducksmiles_compat::ExecuteWithNulls<string_t, int32_t>(                    \
        args.data[0], result, args.size(),                                     \
        [](string_t input, ValidityMask &mask, idx_t idx) -> int32_t {         \
          int32_t val =                                                        \
              RustFunc((const uint8_t *)input.GetData(), input.GetSize(), 0);  \
          if (val < 0) {                                                       \
            mask.SetInvalid(idx);                                              \
            return 0;                                                          \
          }                                                                    \
          return val;                                                          \
        });                                                                    \
  }

DEFINE_STRUCTURE_INT_FUNC(StructureAtomCountFunc, ds_structure_atom_count)
DEFINE_STRUCTURE_INT_FUNC(StructureChainCountFunc, ds_structure_chain_count)
DEFINE_STRUCTURE_INT_FUNC(StructureResidueCountFunc, ds_structure_residue_count)
DEFINE_STRUCTURE_INT_FUNC(StructureModelCountFunc, ds_structure_model_count)

#define DEFINE_STRUCTURE_DOUBLE_FUNC(FuncName, RustFunc)                       \
  static void FuncName(DataChunk &args, ExpressionState &state,                \
                       Vector &result) {                                       \
    ducksmiles_compat::ExecuteWithNulls<string_t, double>(                     \
        args.data[0], result, args.size(),                                     \
        [](string_t input, ValidityMask &mask, idx_t idx) -> double {          \
          double val =                                                         \
              RustFunc((const uint8_t *)input.GetData(), input.GetSize(), 0);  \
          if (std::isnan(val)) {                                               \
            mask.SetInvalid(idx);                                              \
            return 0.0;                                                        \
          }                                                                    \
          return val;                                                          \
        });                                                                    \
  }

DEFINE_STRUCTURE_DOUBLE_FUNC(StructureCentroidXFunc, ds_structure_centroid_x)
DEFINE_STRUCTURE_DOUBLE_FUNC(StructureCentroidYFunc, ds_structure_centroid_y)
DEFINE_STRUCTURE_DOUBLE_FUNC(StructureCentroidZFunc, ds_structure_centroid_z)
DEFINE_STRUCTURE_DOUBLE_FUNC(StructureRadiusOfGyrationFunc,
                             ds_structure_radius_of_gyration)
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

// Both supported DuckDB APIs expose SetFallible(). DuckDB 2 enforces this
// metadata for runtime errors, and TRY/filter pushdown rely on it in both.
static ScalarFunction
FallibleScalarFunction(const char *name, const vector<LogicalType> &arguments,
                       const LogicalType &return_type,
                       scalar_function_t callback) {
  ScalarFunction function(name, arguments, return_type, callback);
  function.SetFallible();
  return function;
}

// Keep overload resolution and conflict behavior identical to
// RegisterFunction(ScalarFunction). Parameter types associate documentation
// with its exact overload in duckdb_functions().
static void RegisterDocumentedScalar(ExtensionLoader &loader,
                                     ScalarFunction function,
                                     vector<string> parameter_names,
                                     const string &description,
                                     const string &example,
                                     const string &category) {
  D_ASSERT(parameter_names.size() == function.arguments.size());
  FunctionDescription documentation;
  documentation.parameter_types = function.arguments;
  documentation.parameter_names = std::move(parameter_names);
  documentation.description = description;
  documentation.examples = {example};
  documentation.categories = {"ducksmiles", category};
  CreateScalarFunctionInfo info(std::move(function));
  info.on_conflict = OnCreateConflict::ALTER_ON_CONFLICT;
  info.descriptions.push_back(std::move(documentation));
  loader.RegisterFunction(std::move(info));
}

static void RegisterDucksmilesFunctions(ExtensionLoader &loader) {
  // --- SMILES ---
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_is_valid", {LogicalType::VARCHAR},
                     LogicalType::BOOLEAN, MolIsValidFunc),
      {"smiles"},
      "Test whether the SMILES is accepted by the supported parser.",
      "mol_is_valid('CCO')", "molecular");
  RegisterDocumentedScalar(loader,
                           ScalarFunction("mol_formula", {LogicalType::VARCHAR},
                                          LogicalType::VARCHAR, MolFormulaFunc),
                           {"smiles"},
                           "Return the molecular formula in Hill order, "
                           "including implicit hydrogens.",
                           "mol_formula('CCO')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_num_atoms", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, MolNumAtomsFunc),
      {"smiles"}, "Count non-hydrogen atoms in a SMILES molecule.",
      "mol_num_atoms('CCO')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_num_fragments", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, MolNumFragmentsFunc),
      {"smiles"}, "Count connected components in the molecular graph.",
      "mol_num_fragments('CCO')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_num_bonds", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, MolNumBondsFunc),
      {"smiles"},
      "Count bonds in a SMILES molecule, with each bond counted once "
      "regardless of order.",
      "mol_num_bonds('CCO')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_formal_charge", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, MolFormalChargeFunc),
      {"smiles"}, "Sum formal charges on all atoms in the molecular graph.",
      "mol_formal_charge('CCO')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_num_explicit_h", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, MolNumExplicitHFunc),
      {"smiles"},
      "Count hydrogen vertices and explicit hydrogen annotations on "
      "non-hydrogen atoms.",
      "mol_num_explicit_h('CCO')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_num_implicit_h", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, MolNumImplicitHFunc),
      {"smiles"},
      "Count implicit hydrogens on non-hydrogen atoms after native aromaticity "
      "perception.",
      "mol_num_implicit_h('CCO')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_num_total_h", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, MolNumTotalHFunc),
      {"smiles"},
      "Count explicit and implicit hydrogens in the molecular graph.",
      "mol_num_total_h('CCO')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_num_single_bonds", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, MolNumSingleBondsFunc),
      {"smiles"},
      "Count single bonds, including bonds to explicit hydrogen atoms.",
      "mol_num_single_bonds('CCO')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_num_double_bonds", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, MolNumDoubleBondsFunc),
      {"smiles"}, "Count double bonds in the molecular graph.",
      "mol_num_double_bonds('CCO')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_num_triple_bonds", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, MolNumTripleBondsFunc),
      {"smiles"}, "Count triple bonds in the molecular graph.",
      "mol_num_triple_bonds('CCO')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_num_aromatic_bonds", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, MolNumAromaticBondsFunc),
      {"smiles"}, "Count bonds classified as aromatic by native perception.",
      "mol_num_aromatic_bonds('c1ccccc1O')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_num_ring_atoms", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, MolNumRingAtomsFunc),
      {"smiles"}, "Count atoms belonging to at least one perceived ring.",
      "mol_num_ring_atoms('c1ccccc1O')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_num_ring_bonds", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, MolNumRingBondsFunc),
      {"smiles"}, "Count bonds belonging to at least one perceived ring.",
      "mol_num_ring_bonds('c1ccccc1O')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_largest_ring_size", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, MolLargestRingSizeFunc),
      {"smiles"},
      "Return the largest ring size in the perceived symmetrized SSSR, or zero "
      "for acyclic molecules.",
      "mol_largest_ring_size('c1ccccc1O')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_num_aromatic_atoms", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, MolNumAromaticAtomsFunc),
      {"smiles"}, "Count atoms classified as aromatic by native perception.",
      "mol_num_aromatic_atoms('c1ccccc1O')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_num_carbons", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, MolNumCarbonsFunc),
      {"smiles"}, "Count carbon atoms in the molecular graph.",
      "mol_num_carbons('CCO')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_num_nitrogens", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, MolNumNitrogensFunc),
      {"smiles"}, "Count nitrogen atoms in the molecular graph.",
      "mol_num_nitrogens('CCO')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_num_oxygens", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, MolNumOxygensFunc),
      {"smiles"}, "Count oxygen atoms in the molecular graph.",
      "mol_num_oxygens('CCO')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_num_halogens", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, MolNumHalogensFunc),
      {"smiles"}, "Count fluorine, chlorine, bromine and iodine atoms.",
      "mol_num_halogens('CCO')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_heteroatom_fraction", {LogicalType::VARCHAR},
                     LogicalType::DOUBLE, MolHeteroatomFractionFunc),
      {"smiles"},
      "Return the fraction of heavy atoms that are neither carbon nor "
      "hydrogen; zero when no heavy atoms exist.",
      "mol_heteroatom_fraction('CCO')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_aromatic_fraction", {LogicalType::VARCHAR},
                     LogicalType::DOUBLE, MolAromaticFractionFunc),
      {"smiles"},
      "Return the fraction of heavy atoms classified as aromatic; zero when no "
      "heavy atoms exist.",
      "mol_aromatic_fraction('c1ccccc1O')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_heavy_atom_mass", {LogicalType::VARCHAR},
                     LogicalType::DOUBLE, MolHeavyAtomMassFunc),
      {"smiles"},
      "Sum isotope-aware heavy-atom masses, using average weights for "
      "unlabelled atoms and excluding all hydrogen isotopes.",
      "mol_heavy_atom_mass('CCO')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_mean_degree", {LogicalType::VARCHAR},
                     LogicalType::DOUBLE, MolMeanDegreeFunc),
      {"smiles"},
      "Return twice the heavy-heavy bond count divided by the heavy-atom "
      "count; zero when no heavy atoms exist.",
      "mol_mean_degree('CCO')", "molecular");
  RegisterDocumentedScalar(loader,
                           ScalarFunction("mol_weight", {LogicalType::VARCHAR},
                                          LogicalType::DOUBLE, MolWeightFunc),
                           {"smiles"},
                           "Calculate molecular weight from standard atomic "
                           "weights, including implicit hydrogens.",
                           "mol_weight('CCO')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_exact_mass", {LogicalType::VARCHAR},
                     LogicalType::DOUBLE, MolExactMassFunc),
      {"smiles"}, "Calculate monoisotopic molecular mass from SMILES.",
      "mol_exact_mass('CCO')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("logp_crippen", {LogicalType::VARCHAR},
                     LogicalType::DOUBLE, LogpCrippenFunc),
      {"smiles"},
      "Estimate octanol/water logP using Wildman-Crippen atom contributions.",
      "logp_crippen('CCO')", "molecular");
  RegisterDocumentedScalar(loader,
                           ScalarFunction("tpsa", {LogicalType::VARCHAR},
                                          LogicalType::DOUBLE, TpsaFunc),
                           {"smiles"},
                           "Calculate topological polar surface area from "
                           "nitrogen and oxygen atom contributions.",
                           "tpsa('CCO')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("canonical_smiles", {LogicalType::VARCHAR},
                     LogicalType::VARCHAR, CanonicalSmilesFunc),
      {"smiles"},
      "Return a deterministic normalized SMILES for the supported molecular "
      "graph subset.",
      "canonical_smiles('CCO')", "molecular");
  RegisterDocumentedScalar(
      loader,
      FallibleScalarFunction("murcko_scaffold", {LogicalType::VARCHAR},
                             LogicalType::VARCHAR, MurckoScaffoldFunc),
      {"smiles"},
      "Extract the Bemis-Murcko ring-and-linker scaffold as SMILES.",
      "murcko_scaffold('Cc1ccccc1')", "molecular");
  RegisterDocumentedScalar(
      loader,
      FallibleScalarFunction("generic_scaffold", {LogicalType::VARCHAR},
                             LogicalType::VARCHAR, GenericScaffoldFunc),
      {"smiles"},
      "Extract a scaffold with carbon atoms and single bonds as SMILES.",
      "generic_scaffold('Cc1ccccc1')", "molecular");
  RegisterDocumentedScalar(
      loader,
      FallibleScalarFunction("ring_systems_json", {LogicalType::VARCHAR},
                             LogicalType::VARCHAR, RingSystemsJsonFunc),
      {"smiles"},
      "Return ring systems and their 1-based atom and bond indices as JSON.",
      "ring_systems_json('Cc1ccccc1')", "molecular");
  RegisterDocumentedScalar(
      loader,
      FallibleScalarFunction("mol_hash",
                             {LogicalType::VARCHAR, LogicalType::VARCHAR},
                             LogicalType::VARCHAR, MolHashFunc),
      {"smiles", "method"},
      "Return a molecular grouping hash for a supported method; list methods "
      "with mol_hash_methods().",
      "mol_hash('CCO', 'element_graph')", "molecular");
  RegisterDocumentedScalar(
      loader,
      FallibleScalarFunction("mol_hash_methods", {}, LogicalType::VARCHAR,
                             MolHashMethodsJsonFunc),
      {}, "List supported molecular hash method names as JSON.",
      "mol_hash_methods()", "molecular");
  RegisterDocumentedScalar(
      loader,
      FallibleScalarFunction("largest_fragment", {LogicalType::VARCHAR},
                             LogicalType::VARCHAR, LargestFragmentFunc),
      {"smiles"}, "Keep the largest connected fragment and return its SMILES.",
      "largest_fragment('CCO')", "molecular");
  RegisterDocumentedScalar(
      loader,
      FallibleScalarFunction("strip_salts", {LogicalType::VARCHAR},
                             LogicalType::VARCHAR, StripSaltsFunc),
      {"smiles"},
      "Remove recognized salt fragments and return the remaining SMILES.",
      "strip_salts('CCO')", "molecular");
  RegisterDocumentedScalar(
      loader,
      FallibleScalarFunction("neutralize_charges", {LogicalType::VARCHAR},
                             LogicalType::VARCHAR, NeutralizeChargesFunc),
      {"smiles"},
      "Neutralize supported charged atom patterns and return SMILES.",
      "neutralize_charges('CCO')", "molecular");
  RegisterDocumentedScalar(
      loader,
      FallibleScalarFunction("normalize_smiles", {LogicalType::VARCHAR},
                             LogicalType::VARCHAR, NormalizeSmilesFunc),
      {"smiles"},
      "Normalize supported functional-group patterns and return SMILES.",
      "normalize_smiles('CCO')", "molecular");
  RegisterDocumentedScalar(
      loader,
      FallibleScalarFunction("fragment_parent", {LogicalType::VARCHAR},
                             LogicalType::VARCHAR, FragmentParentFunc),
      {"smiles"},
      "Normalize, remove salts, select the largest fragment and neutralize "
      "supported charges.",
      "fragment_parent('CCO')", "molecular");
  RegisterDocumentedScalar(
      loader,
      FallibleScalarFunction("mcs_smarts",
                             {LogicalType::VARCHAR, LogicalType::VARCHAR},
                             LogicalType::VARCHAR, McsSmartsFunc),
      {"smiles_a", "smiles_b"},
      "Return a bounded maximum common substructure search result as SMARTS.",
      "mcs_smarts('CCO', 'CCN')", "substructure");
  RegisterDocumentedScalar(
      loader,
      FallibleScalarFunction("mcs_json",
                             {LogicalType::VARCHAR, LogicalType::VARCHAR},
                             LogicalType::VARCHAR, McsJsonFunc),
      {"smiles_a", "smiles_b"},
      "Return a bounded maximum common substructure result and atom mappings "
      "as JSON.",
      "mcs_json('CCO', 'CCN')", "substructure");
  RegisterDocumentedScalar(
      loader,
      FallibleScalarFunction("scaffold_network_json", {LogicalType::VARCHAR},
                             LogicalType::VARCHAR, ScaffoldNetworkJsonFunc),
      {"smiles"}, "Return a bounded ring-removal scaffold network as JSON.",
      "scaffold_network_json('Cc1ccccc1')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("num_h_acceptors", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, NumHAcceptorsFunc),
      {"smiles"}, "Count hydrogen-bond acceptors using the native atom rules.",
      "num_h_acceptors('CCO')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("num_h_donors", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, NumHDonorsFunc),
      {"smiles"},
      "Count hydrogen-bond donor atoms using the native atom rules.",
      "num_h_donors('CCO')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("num_rotatable_bonds", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, NumRotatableBondsFunc),
      {"smiles"}, "Count rotatable bonds using the native strict bond rules.",
      "num_rotatable_bonds('CCO')", "molecular");
  RegisterDocumentedScalar(loader,
                           ScalarFunction("ring_count", {LogicalType::VARCHAR},
                                          LogicalType::INTEGER, RingCountFunc),
                           {"smiles"},
                           "Count rings in the native perceived ring set.",
                           "ring_count('Cc1ccccc1')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("num_aromatic_rings", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, NumAromaticRingsFunc),
      {"smiles"}, "Count fully aromatic rings in the perceived ring set.",
      "num_aromatic_rings('Cc1ccccc1')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("num_heteroatoms", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, NumHeteroatomsFunc),
      {"smiles"}, "Count atoms other than carbon and hydrogen.",
      "num_heteroatoms('CCO')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("fraction_csp3", {LogicalType::VARCHAR},
                     LogicalType::DOUBLE, FractionCsp3Func),
      {"smiles"}, "Return the fraction of carbon atoms classified as sp3.",
      "fraction_csp3('CCO')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_mr", {LogicalType::VARCHAR}, LogicalType::DOUBLE,
                     MolMrFunc),
      {"smiles"},
      "Estimate molar refractivity using Wildman-Crippen atom contributions.",
      "mol_mr('CCO')", "molecular");
  RegisterDocumentedScalar(loader,
                           ScalarFunction("qed", {LogicalType::VARCHAR},
                                          LogicalType::DOUBLE, QedFunc),
                           {"smiles"},
                           "Calculate the weighted QED drug-likeness score "
                           "using native descriptors.",
                           "qed('CCO')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("num_aliphatic_rings", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, NumAliphaticRingsFunc),
      {"smiles"}, "Count rings that are not fully aromatic.",
      "num_aliphatic_rings('Cc1ccccc1')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("num_saturated_rings", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, NumSaturatedRingsFunc),
      {"smiles"}, "Count rings containing only single bonds.",
      "num_saturated_rings('Cc1ccccc1')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("num_aromatic_heterocycles", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, NumAromaticHeterocyclesFunc),
      {"smiles"}, "Count aromatic rings containing a non-carbon atom.",
      "num_aromatic_heterocycles('Cc1ccccc1')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("num_aromatic_carbocycles", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, NumAromaticCarbocyclesFunc),
      {"smiles"}, "Count aromatic rings containing only carbon atoms.",
      "num_aromatic_carbocycles('Cc1ccccc1')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("num_saturated_heterocycles", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, NumSaturatedHeterocyclesFunc),
      {"smiles"}, "Count single-bond rings containing a non-carbon atom.",
      "num_saturated_heterocycles('Cc1ccccc1')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("num_saturated_carbocycles", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, NumSaturatedCarbocyclesFunc),
      {"smiles"}, "Count single-bond rings containing only carbon atoms.",
      "num_saturated_carbocycles('Cc1ccccc1')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("num_aliphatic_heterocycles", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, NumAliphaticHeterocyclesFunc),
      {"smiles"}, "Count non-aromatic rings containing a non-carbon atom.",
      "num_aliphatic_heterocycles('Cc1ccccc1')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("num_aliphatic_carbocycles", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, NumAliphaticCarbocyclesFunc),
      {"smiles"}, "Count non-aromatic rings containing only carbon atoms.",
      "num_aliphatic_carbocycles('Cc1ccccc1')", "molecular");
  // ADMET / drug-likeness rule panels + toxicophore structural alerts
  RegisterDocumentedScalar(
      loader,
      FallibleScalarFunction("admet_json", {LogicalType::VARCHAR},
                             LogicalType::VARCHAR, AdmetJsonFunc),
      {"smiles"},
      "Return native physicochemical descriptors, drug-likeness rule panels "
      "and structural alerts as JSON.",
      "admet_json('CCO')", "molecular");
  RegisterDocumentedScalar(
      loader,
      FallibleScalarFunction("structural_alerts_json", {LogicalType::VARCHAR},
                             LogicalType::VARCHAR, StructuralAlertsJsonFunc),
      {"smiles"}, "Return matched toxicophore structural alerts as JSON.",
      "structural_alerts_json('CCO')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("structural_alert_count", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, StructuralAlertCountFunc),
      {"smiles"}, "Count matched toxicophore structural alerts.",
      "structural_alert_count('CCO')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("lipinski_violations", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, LipinskiViolationsFunc),
      {"smiles"}, "Count violations of the Lipinski rule-of-five thresholds.",
      "lipinski_violations('CCO')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("druglikeness_pass",
                     {LogicalType::VARCHAR, LogicalType::VARCHAR},
                     LogicalType::INTEGER, DruglikenessPassFunc),
      {"smiles", "rule"},
      "Return 1 or 0 for a named drug-likeness rule panel, or NULL for an "
      "unsupported rule or input.",
      "druglikeness_pass('CCO', 'lipinski')", "molecular");
  RegisterDocumentedScalar(
      loader,
      FallibleScalarFunction("pdb_to_pdbqt", {LogicalType::VARCHAR},
                             LogicalType::VARCHAR, PdbToPdbqtFunc),
      {"pdb"}, "Convert protein PDB text to PDBQT with native atom typing.",
      "pdb_to_pdbqt('ATOM      1  CA  ALA A   1       0.000   0.000   0.000  "
      "1.00 20.00           C')",
      "docking");
  // Docking pipeline
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("smiles_to_pdbqt",
                     {LogicalType::VARCHAR, LogicalType::BIGINT},
                     LogicalType::VARCHAR, SmilesToPdbqtFunc),
      {"smiles", "seed"},
      "Generate a seeded 3D ligand conformer from SMILES and return PDBQT.",
      "smiles_to_pdbqt('CCO', 42)", "docking");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("dock",
                     {LogicalType::VARCHAR, LogicalType::VARCHAR,
                      LogicalType::DOUBLE, LogicalType::DOUBLE,
                      LogicalType::DOUBLE, LogicalType::DOUBLE,
                      LogicalType::DOUBLE, LogicalType::DOUBLE,
                      LogicalType::INTEGER, LogicalType::BIGINT},
                     LogicalType::VARCHAR, DockFunc),
      {"smiles", "pdb", "center_x", "center_y", "center_z", "half_extent_x",
       "half_extent_y", "half_extent_z", "n_runs", "seed"},
      "Dock a flexible ligand against a PDB receptor at pH 7.4; return scored "
      "poses as JSON using box half-extents "
      "in angstroms.",
      "dock('CCO', 'ATOM      1  CA  ALA A   1       0.000   0.000   0.000  "
      "1.00 20.00           C', 0.0, 0.0, 0.0, "
      "3.0, 3.0, 3.0, 1, 42)",
      "docking");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction(
          "dock",
          {LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::DOUBLE,
           LogicalType::DOUBLE, LogicalType::DOUBLE, LogicalType::DOUBLE,
           LogicalType::DOUBLE, LogicalType::DOUBLE, LogicalType::INTEGER,
           LogicalType::BIGINT, LogicalType::DOUBLE},
          LogicalType::VARCHAR, DockFunc),
      {"smiles", "pdb", "center_x", "center_y", "center_z", "half_extent_x",
       "half_extent_y", "half_extent_z", "n_runs", "seed", "ph"},
      "Dock a flexible ligand against a PDB receptor at the requested pH; "
      "return scored poses as JSON using box "
      "half-extents in angstroms.",
      "dock('CCO', 'ATOM      1  CA  ALA A   1       0.000   0.000   0.000  "
      "1.00 20.00           C', 0.0, 0.0, 0.0, "
      "3.0, 3.0, 3.0, 1, 42, 7.4)",
      "docking");
  RegisterDocumentedScalar(
      loader,
      FallibleScalarFunction("prepare_receptor",
                             {LogicalType::VARCHAR, LogicalType::DOUBLE},
                             LogicalType::VARCHAR, PrepareReceptorFunc),
      {"pdb", "ph"},
      "Prepare a PDB receptor with native pH-dependent protonation and polar "
      "hydrogens, returning PDBQT.",
      "prepare_receptor('ATOM      1  CA  ALA A   1       0.000   0.000   "
      "0.000  1.00 20.00           C', 7.4)",
      "docking");
  // Virtual-screening benchmark metrics (LIST<DOUBLE> scores, LIST<BOOLEAN>
  // labels)
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("roc_auc",
                     {LogicalType::LIST(LogicalType::DOUBLE),
                      LogicalType::LIST(LogicalType::BOOLEAN)},
                     LogicalType::DOUBLE, RocAucFunc),
      {"scores", "labels"},
      "Calculate ROC AUC for screening scores where lower is better and true "
      "labels mark actives.",
      "roc_auc([-3.0, -2.0, -1.0, 0.0], [true, false, true, false])",
      "benchmark");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("enrichment_factor",
                     {LogicalType::LIST(LogicalType::DOUBLE),
                      LogicalType::LIST(LogicalType::BOOLEAN),
                      LogicalType::DOUBLE},
                     LogicalType::DOUBLE, EnrichmentFactorFunc),
      {"scores", "labels", "fraction"},
      "Calculate active enrichment in the top fraction of a screen ranked by "
      "ascending score.",
      "enrichment_factor([-3.0, -2.0, -1.0, 0.0], [true, false, true, false], "
      "0.5)",
      "benchmark");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("bedroc",
                     {LogicalType::LIST(LogicalType::DOUBLE),
                      LogicalType::LIST(LogicalType::BOOLEAN),
                      LogicalType::DOUBLE},
                     LogicalType::DOUBLE, BedrocFunc),
      {"scores", "labels", "alpha"},
      "Calculate BEDROC early-recognition performance for ascending screening "
      "scores and true active labels.",
      "bedroc([-3.0, -2.0, -1.0, 0.0], [true, false, true, false], 20.0)",
      "benchmark");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_has_substructure",
                     {LogicalType::VARCHAR, LogicalType::VARCHAR},
                     LogicalType::BOOLEAN, MolHasSubstructureFunc),
      {"smiles", "smarts"},
      "Test whether a molecule matches a supported SMARTS pattern.",
      "mol_has_substructure('CCO', 'CO')", "substructure");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_substructure_count",
                     {LogicalType::VARCHAR, LogicalType::VARCHAR},
                     LogicalType::INTEGER, MolSubstructureCountFunc),
      {"smiles", "smarts"},
      "Count unique atom-set matches for a supported SMARTS pattern.",
      "mol_substructure_count('CCO', 'CO')", "substructure");
  RegisterDocumentedScalar(
      loader,
      FallibleScalarFunction("mol_substructure_matches_json",
                             {LogicalType::VARCHAR, LogicalType::VARCHAR},
                             LogicalType::VARCHAR,
                             MolSubstructureMatchesJsonFunc),
      {"smiles", "smarts"},
      "Return unique SMARTS atom-set matches as 1-based atom indices in JSON.",
      "mol_substructure_matches_json('CCO', 'CO')", "substructure");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("add_hydrogens", {LogicalType::VARCHAR},
                     LogicalType::VARCHAR, AddHydrogensFunc),
      {"smiles"},
      "Return SMILES with implicit hydrogens expanded to explicit atoms.",
      "add_hydrogens('CCO')", "molecular");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("morgan_fp_bits", {LogicalType::VARCHAR},
                     LogicalType::BLOB, MorganFpBitsFunc1),
      {"smiles"},
      "Return a native Morgan bit fingerprint with radius 2 and 2048 bits.",
      "morgan_fp_bits('CCO')", "fingerprint");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction(
          "morgan_fp_bits",
          {LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::INTEGER},
          LogicalType::BLOB, MorganFpBitsFunc3),
      {"smiles", "radius", "n_bits"},
      "Return a native Morgan bit fingerprint with the requested radius and "
      "bit count.",
      "morgan_fp_bits('CCO', 2, 2048)", "fingerprint");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("maccs_keys", {LogicalType::VARCHAR}, LogicalType::BLOB,
                     MaccsKeysFunc),
      {"smiles"},
      "Return the 166 MACCS structural keys in a 21-byte bit fingerprint.",
      "maccs_keys('CCO')", "fingerprint");
  RegisterDocumentedScalar(
      loader,
      FallibleScalarFunction("tanimoto_bit",
                             {LogicalType::BLOB, LogicalType::BLOB},
                             LogicalType::DOUBLE, TanimotoBitFunc),
      {"fingerprint_a", "fingerprint_b"},
      "Calculate Tanimoto similarity between equal-length BLOB bit "
      "fingerprints.",
      "tanimoto_bit(morgan_fp_bits('CCO'), morgan_fp_bits('CCN'))",
      "similarity");
  RegisterDocumentedScalar(
      loader,
      FallibleScalarFunction("dice_bit", {LogicalType::BLOB, LogicalType::BLOB},
                             LogicalType::DOUBLE, DiceBitFunc),
      {"fingerprint_a", "fingerprint_b"},
      "Calculate Dice similarity between equal-length BLOB bit fingerprints.",
      "dice_bit(morgan_fp_bits('CCO'), morgan_fp_bits('CCN'))", "similarity");
  RegisterDocumentedScalar(
      loader,
      FallibleScalarFunction("cosine_bit",
                             {LogicalType::BLOB, LogicalType::BLOB},
                             LogicalType::DOUBLE, CosineBitFunc),
      {"fingerprint_a", "fingerprint_b"},
      "Calculate cosine similarity between equal-length BLOB bit fingerprints.",
      "cosine_bit(morgan_fp_bits('CCO'), morgan_fp_bits('CCN'))", "similarity");
  RegisterDocumentedScalar(
      loader,
      FallibleScalarFunction("kulczynski_bit",
                             {LogicalType::BLOB, LogicalType::BLOB},
                             LogicalType::DOUBLE, KulczynskiBitFunc),
      {"fingerprint_a", "fingerprint_b"},
      "Calculate Kulczynski similarity between equal-length BLOB bit "
      "fingerprints.",
      "kulczynski_bit(morgan_fp_bits('CCO'), morgan_fp_bits('CCN'))",
      "similarity");
  RegisterDocumentedScalar(
      loader,
      FallibleScalarFunction("sokal_bit",
                             {LogicalType::BLOB, LogicalType::BLOB},
                             LogicalType::DOUBLE, SokalBitFunc),
      {"fingerprint_a", "fingerprint_b"},
      "Calculate Sokal similarity between equal-length BLOB bit fingerprints.",
      "sokal_bit(morgan_fp_bits('CCO'), morgan_fp_bits('CCN'))", "similarity");
  RegisterDocumentedScalar(
      loader,
      FallibleScalarFunction("mcconnaughey_bit",
                             {LogicalType::BLOB, LogicalType::BLOB},
                             LogicalType::DOUBLE, McConnaugheyBitFunc),
      {"fingerprint_a", "fingerprint_b"},
      "Calculate McConnaughey similarity between equal-length BLOB bit "
      "fingerprints.",
      "mcconnaughey_bit(morgan_fp_bits('CCO'), morgan_fp_bits('CCN'))",
      "similarity");
  RegisterDocumentedScalar(
      loader,
      FallibleScalarFunction("asymmetric_bit",
                             {LogicalType::BLOB, LogicalType::BLOB},
                             LogicalType::DOUBLE, AsymmetricBitFunc),
      {"fingerprint_a", "fingerprint_b"},
      "Calculate asymmetric (intersection divided by the smaller on-bit count) "
      "similarity "
      "between equal-length BLOB bit fingerprints.",
      "asymmetric_bit(morgan_fp_bits('CCO'), morgan_fp_bits('CCN'))",
      "similarity");
  RegisterDocumentedScalar(
      loader,
      FallibleScalarFunction("braun_blanquet_bit",
                             {LogicalType::BLOB, LogicalType::BLOB},
                             LogicalType::DOUBLE, BraunBlanquetBitFunc),
      {"fingerprint_a", "fingerprint_b"},
      "Calculate Braun-Blanquet similarity between equal-length BLOB bit "
      "fingerprints.",
      "braun_blanquet_bit(morgan_fp_bits('CCO'), morgan_fp_bits('CCN'))",
      "similarity");
  RegisterDocumentedScalar(
      loader,
      FallibleScalarFunction("russel_bit",
                             {LogicalType::BLOB, LogicalType::BLOB},
                             LogicalType::DOUBLE, RusselBitFunc),
      {"fingerprint_a", "fingerprint_b"},
      "Calculate Russel similarity between equal-length BLOB bit fingerprints.",
      "russel_bit(morgan_fp_bits('CCO'), morgan_fp_bits('CCN'))", "similarity");
  RegisterDocumentedScalar(
      loader,
      FallibleScalarFunction("tversky_bit",
                             {LogicalType::BLOB, LogicalType::BLOB,
                              LogicalType::DOUBLE, LogicalType::DOUBLE},
                             LogicalType::DOUBLE, TverskyBitFunc),
      {"fingerprint_a", "fingerprint_b", "alpha", "beta"},
      "Calculate Tversky similarity between equal-length BLOB bit fingerprints "
      "with alpha and beta weights.",
      "tversky_bit(morgan_fp_bits('CCO'), morgan_fp_bits('CCN'), 0.5, 0.5)",
      "similarity");

  // --- InChI layer extraction ---
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("inchi_is_valid", {LogicalType::VARCHAR},
                     LogicalType::BOOLEAN, InchiIsValidFunc),
      {"inchi"}, "Test whether text has a supported InChI syntax.",
      "inchi_is_valid('InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)')", "inchi");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("inchi_is_standard", {LogicalType::VARCHAR},
                     LogicalType::BOOLEAN, InchiIsStandardFunc),
      {"inchi"}, "Test whether an InChI has the standard 1S version marker.",
      "inchi_is_standard('InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)')", "inchi");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("inchi_version", {LogicalType::VARCHAR},
                     LogicalType::VARCHAR, InchiVersionFunc),
      {"inchi"},
      "Extract the InChI version marker, including S for standard InChI.",
      "inchi_version('InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)')", "inchi");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("inchi_formula", {LogicalType::VARCHAR},
                     LogicalType::VARCHAR, InchiFormulaFunc),
      {"inchi"}, "Extract the molecular formula layer from an InChI.",
      "inchi_formula('InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)')", "inchi");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("inchi_connections", {LogicalType::VARCHAR},
                     LogicalType::VARCHAR, InchiConnectionsFunc),
      {"inchi"},
      "Extract the InChI connectivity layer (c), or an empty string when "
      "absent.",
      "inchi_connections('InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)')", "inchi");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("inchi_hydrogens", {LogicalType::VARCHAR},
                     LogicalType::VARCHAR, InchiHydrogensFunc),
      {"inchi"},
      "Extract the InChI hydrogen layer (h), or an empty string when absent.",
      "inchi_hydrogens('InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)')", "inchi");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("inchi_charge", {LogicalType::VARCHAR},
                     LogicalType::VARCHAR, InchiChargeFunc),
      {"inchi"},
      "Extract the InChI charge layer (q), or an empty string when absent.",
      "inchi_charge('InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)')", "inchi");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("inchi_stereo_bond", {LogicalType::VARCHAR},
                     LogicalType::VARCHAR, InchiStereoBondFunc),
      {"inchi"},
      "Extract the InChI double-bond stereo layer (b), or an empty string when "
      "absent.",
      "inchi_stereo_bond('InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)')", "inchi");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("inchi_stereo_tetrahedral", {LogicalType::VARCHAR},
                     LogicalType::VARCHAR, InchiStereoTetrahedralFunc),
      {"inchi"},
      "Extract the InChI tetrahedral stereo layer (t), or an empty string when "
      "absent.",
      "inchi_stereo_tetrahedral('InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)')",
      "inchi");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("inchi_has_stereo", {LogicalType::VARCHAR},
                     LogicalType::BOOLEAN, InchiHasStereoFunc),
      {"inchi"},
      "Test whether an InChI contains bond or tetrahedral stereo layers.",
      "inchi_has_stereo('InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)')", "inchi");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("inchi_num_stereo_centers", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, InchiNumStereoCentersFunc),
      {"inchi"}, "Count entries in the InChI tetrahedral stereo layer.",
      "inchi_num_stereo_centers('InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)')",
      "inchi");
  // --- InChIKey ---
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("inchikey_is_valid", {LogicalType::VARCHAR},
                     LogicalType::BOOLEAN, InchikeyIsValidFunc),
      {"inchikey"},
      "Test the uppercase 14-10-1 character format of an InChIKey.",
      "inchikey_is_valid('QTBSBXVTEAMEQO-UHFFFAOYSA-N')", "inchi");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("inchikey_connectivity", {LogicalType::VARCHAR},
                     LogicalType::VARCHAR, InchikeyConnectivityFunc),
      {"inchikey"},
      "Extract the first 14-character connectivity block of an InChIKey.",
      "inchikey_connectivity('QTBSBXVTEAMEQO-UHFFFAOYSA-N')", "inchi");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("inchikey_stereo", {LogicalType::VARCHAR},
                     LogicalType::VARCHAR, InchikeyStereofunc),
      {"inchikey"}, "Extract the second 10-character block of an InChIKey.",
      "inchikey_stereo('QTBSBXVTEAMEQO-UHFFFAOYSA-N')", "inchi");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("inchikey_protonation", {LogicalType::VARCHAR},
                     LogicalType::VARCHAR, InchikeyProtonationFunc),
      {"inchikey"}, "Extract the final protonation character of an InChIKey.",
      "inchikey_protonation('QTBSBXVTEAMEQO-UHFFFAOYSA-N')", "inchi");

  // --- Comparison ---
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("inchi_skeleton_match",
                     {LogicalType::VARCHAR, LogicalType::VARCHAR},
                     LogicalType::BOOLEAN, InchiSkeletonMatchFunc),
      {"inchi_a", "inchi_b"},
      "Compare InChI formula, connectivity and hydrogen layers while ignoring "
      "stereo layers.",
      "inchi_skeleton_match('InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)', "
      "'InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)')",
      "inchi");

  // --- MOL/SDF ---
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_block_formula", {LogicalType::VARCHAR},
                     LogicalType::VARCHAR, MolBlockFormulaFunc),
      {"mol_block"},
      "Return the formula from the atoms recorded in the first MOL block.",
      "mol_block_formula('water\n  ducksmiles\n\n  1  0  0  0  0  0  0  0  0  "
      "0999 V2000\n    0.0000    0.0000    "
      "0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n> "
      "<ID>\nwater\n\n$$$$\n')",
      "mol_sdf");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_block_weight", {LogicalType::VARCHAR},
                     LogicalType::DOUBLE, MolBlockWeightFunc),
      {"mol_block"},
      "Calculate molecular weight from the atoms recorded in the first MOL "
      "block.",
      "mol_block_weight('water\n  ducksmiles\n\n  1  0  0  0  0  0  0  0  0  "
      "0999 V2000\n    0.0000    0.0000    "
      "0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n> "
      "<ID>\nwater\n\n$$$$\n')",
      "mol_sdf");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_block_num_atoms", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, MolBlockNumAtomsFunc),
      {"mol_block"}, "Count atoms recorded in the first MOL block.",
      "mol_block_num_atoms('water\n  ducksmiles\n\n  1  0  0  0  0  0  0  0  0 "
      " 0999 V2000\n    0.0000    0.0000    "
      "0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n> "
      "<ID>\nwater\n\n$$$$\n')",
      "mol_sdf");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_block_num_bonds", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, MolBlockNumBondsFunc),
      {"mol_block"}, "Count bonds recorded in the first MOL block.",
      "mol_block_num_bonds('water\n  ducksmiles\n\n  1  0  0  0  0  0  0  0  0 "
      " 0999 V2000\n    0.0000    0.0000    "
      "0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n> "
      "<ID>\nwater\n\n$$$$\n')",
      "mol_sdf");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_block_name", {LogicalType::VARCHAR},
                     LogicalType::VARCHAR, MolBlockNameFunc),
      {"mol_block"}, "Extract the name of the first MOL block.",
      "mol_block_name('water\n  ducksmiles\n\n  1  0  0  0  0  0  0  0  0  "
      "0999 V2000\n    0.0000    0.0000    "
      "0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n> "
      "<ID>\nwater\n\n$$$$\n')",
      "mol_sdf");
  RegisterDocumentedScalar(
      loader,
      FallibleScalarFunction("mol_block_property",
                             {LogicalType::VARCHAR, LogicalType::VARCHAR},
                             LogicalType::VARCHAR, MolBlockPropertyFunc),
      {"mol_block", "property_name"},
      "Extract a named property from the first MOL/SDF record, or NULL when "
      "absent.",
      "mol_block_property('water\n  ducksmiles\n\n  1  0  0  0  0  0  0  0  0  "
      "0999 V2000\n    0.0000    0.0000    "
      "0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n> "
      "<ID>\nwater\n\n$$$$\n', 'ID')",
      "mol_sdf");
  RegisterDocumentedScalar(
      loader,
      FallibleScalarFunction("mol_block_properties_json",
                             {LogicalType::VARCHAR}, LogicalType::VARCHAR,
                             MolBlockPropertiesJsonFunc),
      {"mol_block"}, "Return the first MOL/SDF record properties as JSON.",
      "mol_block_properties_json('water\n  ducksmiles\n\n  1  0  0  0  0  0  0 "
      " 0  0  0999 V2000\n    0.0000    "
      "0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n> "
      "<ID>\nwater\n\n$$$$\n')",
      "mol_sdf");
  RegisterDocumentedScalar(
      loader,
      FallibleScalarFunction("mol_block_atoms_json", {LogicalType::VARCHAR},
                             LogicalType::VARCHAR, MolBlockAtomsJsonFunc),
      {"mol_block"},
      "Return atoms and coordinates from the first MOL block as JSON.",
      "mol_block_atoms_json('water\n  ducksmiles\n\n  1  0  0  0  0  0  0  0  "
      "0  0999 V2000\n    0.0000    0.0000    "
      "0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n> "
      "<ID>\nwater\n\n$$$$\n')",
      "mol_sdf");
  RegisterDocumentedScalar(
      loader,
      FallibleScalarFunction("mol_block_bonds_json", {LogicalType::VARCHAR},
                             LogicalType::VARCHAR, MolBlockBondsJsonFunc),
      {"mol_block"},
      "Return bonds and 1-based atom indices from the first MOL block as JSON.",
      "mol_block_bonds_json('water\n  ducksmiles\n\n  1  0  0  0  0  0  0  0  "
      "0  0999 V2000\n    0.0000    0.0000    "
      "0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n> "
      "<ID>\nwater\n\n$$$$\n')",
      "mol_sdf");
  RegisterDocumentedScalar(
      loader,
      FallibleScalarFunction("mol_block_json", {LogicalType::VARCHAR},
                             LogicalType::VARCHAR, MolBlockJsonFunc),
      {"mol_block"},
      "Return the first MOL record, geometry and properties as JSON.",
      "mol_block_json('water\n  ducksmiles\n\n  1  0  0  0  0  0  0  0  0  "
      "0999 V2000\n    0.0000    0.0000    "
      "0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n> "
      "<ID>\nwater\n\n$$$$\n')",
      "mol_sdf");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_block_has_3d", {LogicalType::VARCHAR},
                     LogicalType::BOOLEAN, MolBlockHas3dFunc),
      {"mol_block"},
      "Test whether any atom in the first MOL block has an absolute z "
      "coordinate above 0.0001.",
      "mol_block_has_3d('water\n  ducksmiles\n\n  1  0  0  0  0  0  0  0  0  "
      "0999 V2000\n    0.0000    0.0000    "
      "0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n> "
      "<ID>\nwater\n\n$$$$\n')",
      "mol_sdf");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_block_centroid_x", {LogicalType::VARCHAR},
                     LogicalType::DOUBLE, MolBlockCentroidXFunc),
      {"mol_block"},
      "Return the mean x coordinate of atoms in the first MOL block.",
      "mol_block_centroid_x('water\n  ducksmiles\n\n  1  0  0  0  0  0  0  0  "
      "0  0999 V2000\n    0.0000    0.0000    "
      "0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n> "
      "<ID>\nwater\n\n$$$$\n')",
      "mol_sdf");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_block_centroid_y", {LogicalType::VARCHAR},
                     LogicalType::DOUBLE, MolBlockCentroidYFunc),
      {"mol_block"},
      "Return the mean y coordinate of atoms in the first MOL block.",
      "mol_block_centroid_y('water\n  ducksmiles\n\n  1  0  0  0  0  0  0  0  "
      "0  0999 V2000\n    0.0000    0.0000    "
      "0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n> "
      "<ID>\nwater\n\n$$$$\n')",
      "mol_sdf");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_block_centroid_z", {LogicalType::VARCHAR},
                     LogicalType::DOUBLE, MolBlockCentroidZFunc),
      {"mol_block"},
      "Return the mean z coordinate of atoms in the first MOL block.",
      "mol_block_centroid_z('water\n  ducksmiles\n\n  1  0  0  0  0  0  0  0  "
      "0  0999 V2000\n    0.0000    0.0000    "
      "0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n> "
      "<ID>\nwater\n\n$$$$\n')",
      "mol_sdf");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_block_radius_of_gyration", {LogicalType::VARCHAR},
                     LogicalType::DOUBLE, MolBlockRadiusOfGyrationFunc),
      {"mol_block"},
      "Calculate the unweighted radius of gyration of coordinates in the first "
      "MOL block.",
      "mol_block_radius_of_gyration('water\n  ducksmiles\n\n  1  0  0  0  0  0 "
      " 0  0  0  0999 V2000\n    0.0000    "
      "0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n> "
      "<ID>\nwater\n\n$$$$\n')",
      "mol_sdf");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_block_min_x", {LogicalType::VARCHAR},
                     LogicalType::DOUBLE, MolBlockMinXFunc),
      {"mol_block"},
      "Return the minimum x coordinate of atoms in the first MOL block.",
      "mol_block_min_x('water\n  ducksmiles\n\n  1  0  0  0  0  0  0  0  0  "
      "0999 V2000\n    0.0000    0.0000    "
      "0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n> "
      "<ID>\nwater\n\n$$$$\n')",
      "mol_sdf");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_block_max_x", {LogicalType::VARCHAR},
                     LogicalType::DOUBLE, MolBlockMaxXFunc),
      {"mol_block"},
      "Return the maximum x coordinate of atoms in the first MOL block.",
      "mol_block_max_x('water\n  ducksmiles\n\n  1  0  0  0  0  0  0  0  0  "
      "0999 V2000\n    0.0000    0.0000    "
      "0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n> "
      "<ID>\nwater\n\n$$$$\n')",
      "mol_sdf");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_block_min_y", {LogicalType::VARCHAR},
                     LogicalType::DOUBLE, MolBlockMinYFunc),
      {"mol_block"},
      "Return the minimum y coordinate of atoms in the first MOL block.",
      "mol_block_min_y('water\n  ducksmiles\n\n  1  0  0  0  0  0  0  0  0  "
      "0999 V2000\n    0.0000    0.0000    "
      "0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n> "
      "<ID>\nwater\n\n$$$$\n')",
      "mol_sdf");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_block_max_y", {LogicalType::VARCHAR},
                     LogicalType::DOUBLE, MolBlockMaxYFunc),
      {"mol_block"},
      "Return the maximum y coordinate of atoms in the first MOL block.",
      "mol_block_max_y('water\n  ducksmiles\n\n  1  0  0  0  0  0  0  0  0  "
      "0999 V2000\n    0.0000    0.0000    "
      "0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n> "
      "<ID>\nwater\n\n$$$$\n')",
      "mol_sdf");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_block_min_z", {LogicalType::VARCHAR},
                     LogicalType::DOUBLE, MolBlockMinZFunc),
      {"mol_block"},
      "Return the minimum z coordinate of atoms in the first MOL block.",
      "mol_block_min_z('water\n  ducksmiles\n\n  1  0  0  0  0  0  0  0  0  "
      "0999 V2000\n    0.0000    0.0000    "
      "0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n> "
      "<ID>\nwater\n\n$$$$\n')",
      "mol_sdf");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("mol_block_max_z", {LogicalType::VARCHAR},
                     LogicalType::DOUBLE, MolBlockMaxZFunc),
      {"mol_block"},
      "Return the maximum z coordinate of atoms in the first MOL block.",
      "mol_block_max_z('water\n  ducksmiles\n\n  1  0  0  0  0  0  0  0  0  "
      "0999 V2000\n    0.0000    0.0000    "
      "0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n> "
      "<ID>\nwater\n\n$$$$\n')",
      "mol_sdf");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("sdf_count", {LogicalType::VARCHAR}, LogicalType::INTEGER,
                     SdfCountFunc),
      {"sdf"}, "Count parsed molecule records in SDF text.",
      "sdf_count('water\n  ducksmiles\n\n  1  0  0  0  0  0  0  0  0  0999 "
      "V2000\n    0.0000    0.0000    0.0000 O   "
      "0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n> <ID>\nwater\n\n$$$$\n')",
      "mol_sdf");
  RegisterDocumentedScalar(
      loader,
      FallibleScalarFunction(
          "sdf_property",
          {LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::VARCHAR},
          LogicalType::VARCHAR, SdfPropertyFunc),
      {"sdf", "record_index", "property_name"},
      "Extract a named SDF property using a 1-based record index, or NULL when "
      "absent.",
      "sdf_property('water\n  ducksmiles\n\n  1  0  0  0  0  0  0  0  0  0999 "
      "V2000\n    0.0000    0.0000    0.0000 "
      "O   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n> "
      "<ID>\nwater\n\n$$$$\n', 1, 'ID')",
      "mol_sdf");
  RegisterDocumentedScalar(
      loader,
      FallibleScalarFunction("sdf_properties_json", {LogicalType::VARCHAR},
                             LogicalType::VARCHAR, SdfPropertiesJsonFunc),
      {"sdf"}, "Return each parsed SDF record name and properties as JSON.",
      "sdf_properties_json('water\n  ducksmiles\n\n  1  0  0  0  0  0  0  0  0 "
      " 0999 V2000\n    0.0000    0.0000    "
      "0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n> "
      "<ID>\nwater\n\n$$$$\n')",
      "mol_sdf");

  // --- PDB/CIF/XYZ structure ---
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("structure_atom_count", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, StructureAtomCountFunc),
      {"structure"}, "Count atoms in auto-detected PDB, mmCIF or XYZ text.",
      "structure_atom_count('ATOM      1  CA  ALA A   1       0.000   0.000   "
      "0.000  1.00 20.00           C')",
      "structure");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("structure_chain_count", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, StructureChainCountFunc),
      {"structure"},
      "Count distinct chains in auto-detected PDB, mmCIF or XYZ text.",
      "structure_chain_count('ATOM      1  CA  ALA A   1       0.000   0.000   "
      "0.000  1.00 20.00           C')",
      "structure");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("structure_residue_count", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, StructureResidueCountFunc),
      {"structure"},
      "Count distinct residues in auto-detected PDB, mmCIF or XYZ text.",
      "structure_residue_count('ATOM      1  CA  ALA A   1       0.000   0.000 "
      "  0.000  1.00 20.00           C')",
      "structure");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("structure_model_count", {LogicalType::VARCHAR},
                     LogicalType::INTEGER, StructureModelCountFunc),
      {"structure"}, "Count models in auto-detected PDB, mmCIF or XYZ text.",
      "structure_model_count('ATOM      1  CA  ALA A   1       0.000   0.000   "
      "0.000  1.00 20.00           C')",
      "structure");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("structure_centroid_x", {LogicalType::VARCHAR},
                     LogicalType::DOUBLE, StructureCentroidXFunc),
      {"structure"},
      "Return the mean x coordinate of atoms in auto-detected PDB, mmCIF or "
      "XYZ text.",
      "structure_centroid_x('ATOM      1  CA  ALA A   1       0.000   0.000   "
      "0.000  1.00 20.00           C')",
      "structure");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("structure_centroid_y", {LogicalType::VARCHAR},
                     LogicalType::DOUBLE, StructureCentroidYFunc),
      {"structure"},
      "Return the mean y coordinate of atoms in auto-detected PDB, mmCIF or "
      "XYZ text.",
      "structure_centroid_y('ATOM      1  CA  ALA A   1       0.000   0.000   "
      "0.000  1.00 20.00           C')",
      "structure");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("structure_centroid_z", {LogicalType::VARCHAR},
                     LogicalType::DOUBLE, StructureCentroidZFunc),
      {"structure"},
      "Return the mean z coordinate of atoms in auto-detected PDB, mmCIF or "
      "XYZ text.",
      "structure_centroid_z('ATOM      1  CA  ALA A   1       0.000   0.000   "
      "0.000  1.00 20.00           C')",
      "structure");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("structure_radius_of_gyration", {LogicalType::VARCHAR},
                     LogicalType::DOUBLE, StructureRadiusOfGyrationFunc),
      {"structure"},
      "Calculate the unweighted radius of gyration of a PDB, mmCIF or XYZ "
      "structure.",
      "structure_radius_of_gyration('ATOM      1  CA  ALA A   1       0.000   "
      "0.000   0.000  "
      "1.00 20.00           C')",
      "structure");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("structure_min_x", {LogicalType::VARCHAR},
                     LogicalType::DOUBLE, StructureMinXFunc),
      {"structure"},
      "Return the minimum x coordinate of atoms in auto-detected PDB, mmCIF or "
      "XYZ text.",
      "structure_min_x('ATOM      1  CA  ALA A   1       0.000   0.000   0.000 "
      " 1.00 20.00           C')",
      "structure");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("structure_max_x", {LogicalType::VARCHAR},
                     LogicalType::DOUBLE, StructureMaxXFunc),
      {"structure"},
      "Return the maximum x coordinate of atoms in auto-detected PDB, mmCIF or "
      "XYZ text.",
      "structure_max_x('ATOM      1  CA  ALA A   1       0.000   0.000   0.000 "
      " 1.00 20.00           C')",
      "structure");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("structure_min_y", {LogicalType::VARCHAR},
                     LogicalType::DOUBLE, StructureMinYFunc),
      {"structure"},
      "Return the minimum y coordinate of atoms in auto-detected PDB, mmCIF or "
      "XYZ text.",
      "structure_min_y('ATOM      1  CA  ALA A   1       0.000   0.000   0.000 "
      " 1.00 20.00           C')",
      "structure");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("structure_max_y", {LogicalType::VARCHAR},
                     LogicalType::DOUBLE, StructureMaxYFunc),
      {"structure"},
      "Return the maximum y coordinate of atoms in auto-detected PDB, mmCIF or "
      "XYZ text.",
      "structure_max_y('ATOM      1  CA  ALA A   1       0.000   0.000   0.000 "
      " 1.00 20.00           C')",
      "structure");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("structure_min_z", {LogicalType::VARCHAR},
                     LogicalType::DOUBLE, StructureMinZFunc),
      {"structure"},
      "Return the minimum z coordinate of atoms in auto-detected PDB, mmCIF or "
      "XYZ text.",
      "structure_min_z('ATOM      1  CA  ALA A   1       0.000   0.000   0.000 "
      " 1.00 20.00           C')",
      "structure");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("structure_max_z", {LogicalType::VARCHAR},
                     LogicalType::DOUBLE, StructureMaxZFunc),
      {"structure"},
      "Return the maximum z coordinate of atoms in auto-detected PDB, mmCIF or "
      "XYZ text.",
      "structure_max_z('ATOM      1  CA  ALA A   1       0.000   0.000   0.000 "
      " 1.00 20.00           C')",
      "structure");

  // --- SELFIES ---
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("smiles_to_selfies", {LogicalType::VARCHAR},
                     LogicalType::VARCHAR, SmilesToSelfiesFunc),
      {"smiles"}, "Encode supported SMILES as SELFIES.",
      "smiles_to_selfies('CCO')", "selfies");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("selfies_to_smiles", {LogicalType::VARCHAR},
                     LogicalType::VARCHAR, SelfiesToSmilesFunc),
      {"selfies"}, "Decode supported SELFIES tokens to SMILES.",
      "selfies_to_smiles('[C][C][O]')", "selfies");
  RegisterDocumentedScalar(
      loader,
      ScalarFunction("selfies_is_valid", {LogicalType::VARCHAR},
                     LogicalType::BOOLEAN, SelfiesIsValidFunc),
      {"selfies"},
      "Test whether a SELFIES string can be decoded by the supported decoder.",
      "selfies_is_valid('[C][C][O]')", "selfies");
}

// ============================================================================
// Extension lifecycle
// ============================================================================

void DucksmilesExtension::Load(ExtensionLoader &loader) {
  RegisterDucksmilesFunctions(loader);
}

std::string DucksmilesExtension::Name() { return "ducksmiles"; }

std::string DucksmilesExtension::Version() const {
#ifdef EXT_VERSION_DUCKSMILES
  return EXT_VERSION_DUCKSMILES;
#else
  return "0.1.0";
#endif
}

} // namespace duckdb

extern "C" {

DUCKDB_EXTENSION_API void
ducksmiles_duckdb_cpp_init(duckdb::ExtensionLoader &loader) {
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
