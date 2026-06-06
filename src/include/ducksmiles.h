#pragma once

#include <cstddef>
#include <cstdint>

// C FFI functions exported from Rust static libraries

extern "C" {

// ===================== SMILES crate =====================

int32_t ds_mol_is_valid(const uint8_t *ptr, size_t len);
int32_t ds_mol_num_atoms(const uint8_t *ptr, size_t len);
int32_t ds_mol_num_bonds(const uint8_t *ptr, size_t len);
int32_t ds_mol_formula(const uint8_t *ptr, size_t len, uint8_t *out, size_t out_cap);
double ds_mol_weight(const uint8_t *ptr, size_t len);
double ds_mol_exact_mass(const uint8_t *ptr, size_t len);
double ds_logp_crippen(const uint8_t *ptr, size_t len);
double ds_tpsa(const uint8_t *ptr, size_t len);

// SMILES with explicit H atoms (verbose bracket form). Returns length written, or -1.
int32_t ds_add_hydrogens(const uint8_t *ptr, size_t len, uint8_t *out, size_t out_cap);

// Morgan/ECFP fingerprint. Writes ceil(n_bits/8) bytes. Returns bytes written or -1.
int32_t ds_morgan_fp_bits(const uint8_t *ptr, size_t len,
                          uint32_t radius, uint32_t n_bits,
                          uint8_t *out, size_t out_cap);

// Tanimoto similarity over two raw fingerprint BLOBs. Returns NaN on length
// mismatch; 0.0 on both-empty; popcount(a & b) / popcount(a | b) otherwise.
double ds_tanimoto_bit(const uint8_t *a_ptr, size_t a_len,
                       const uint8_t *b_ptr, size_t b_len);

// ===================== InChI crate =====================

int32_t ds_inchi_is_valid(const uint8_t *ptr, size_t len);
int32_t ds_inchi_is_standard(const uint8_t *ptr, size_t len);
int32_t ds_inchi_version(const uint8_t *ptr, size_t len, uint8_t *out, size_t out_cap);
int32_t ds_inchi_formula(const uint8_t *ptr, size_t len, uint8_t *out, size_t out_cap);
int32_t ds_inchi_connections(const uint8_t *ptr, size_t len, uint8_t *out, size_t out_cap);
int32_t ds_inchi_hydrogens(const uint8_t *ptr, size_t len, uint8_t *out, size_t out_cap);
int32_t ds_inchi_charge(const uint8_t *ptr, size_t len, uint8_t *out, size_t out_cap);
int32_t ds_inchi_stereo_bond(const uint8_t *ptr, size_t len, uint8_t *out, size_t out_cap);
int32_t ds_inchi_stereo_tetrahedral(const uint8_t *ptr, size_t len, uint8_t *out, size_t out_cap);
int32_t ds_inchi_has_stereo(const uint8_t *ptr, size_t len);
int32_t ds_inchi_num_stereo_centers(const uint8_t *ptr, size_t len);
int32_t ds_inchikey_is_valid(const uint8_t *ptr, size_t len);
int32_t ds_inchikey_connectivity(const uint8_t *ptr, size_t len, uint8_t *out, size_t out_cap);
int32_t ds_inchikey_stereo(const uint8_t *ptr, size_t len, uint8_t *out, size_t out_cap);
int32_t ds_inchikey_protonation(const uint8_t *ptr, size_t len, uint8_t *out, size_t out_cap);
int32_t ds_inchi_skeleton_match(const uint8_t *a_ptr, size_t a_len,
                                const uint8_t *b_ptr, size_t b_len);

// ===================== MOL crate =====================

int32_t ds_sdf_count(const uint8_t *data, size_t len);
int32_t ds_mol_block_formula(const uint8_t *data, size_t len, uint8_t *out, size_t cap);
double ds_mol_block_weight(const uint8_t *data, size_t len);
int32_t ds_mol_block_num_atoms(const uint8_t *data, size_t len);
int32_t ds_mol_block_num_bonds(const uint8_t *data, size_t len);
int32_t ds_mol_block_name(const uint8_t *data, size_t len, uint8_t *out, size_t cap);

// ===================== PDB crate =====================

// format: 0=auto, 1=pdb, 2=cif, 3=xyz
int32_t ds_structure_atom_count(const uint8_t *data, size_t len, uint8_t format);
int32_t ds_structure_chain_count(const uint8_t *data, size_t len, uint8_t format);
int32_t ds_structure_residue_count(const uint8_t *data, size_t len, uint8_t format);
int32_t ds_structure_model_count(const uint8_t *data, size_t len, uint8_t format);

// ===================== SELFIES crate =====================

// SMILES → SELFIES. Returns length written, or -1.
int32_t ds_smiles_to_selfies(const uint8_t *ptr, size_t len, uint8_t *out, size_t cap);

// SELFIES → SMILES. Returns length written, or -1.
int32_t ds_selfies_to_smiles(const uint8_t *ptr, size_t len, uint8_t *out, size_t cap);

// Validate SELFIES: returns 1 if valid, 0 otherwise.
int32_t ds_selfies_is_valid(const uint8_t *ptr, size_t len);

} // extern "C"
