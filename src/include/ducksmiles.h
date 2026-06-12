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
int32_t ds_canonical_smiles(const uint8_t *ptr, size_t len, uint8_t *out, size_t out_cap);
int32_t ds_murcko_scaffold(const uint8_t *ptr, size_t len, uint8_t *out, size_t out_cap);
int32_t ds_generic_scaffold(const uint8_t *ptr, size_t len, uint8_t *out, size_t out_cap);
int32_t ds_ring_systems_json(const uint8_t *ptr, size_t len, uint8_t *out, size_t out_cap);
int32_t ds_mol_hash(const uint8_t *smiles_ptr, size_t smiles_len,
                    const uint8_t *method_ptr, size_t method_len,
                    uint8_t *out, size_t out_cap);
int32_t ds_mol_hash_methods_json(uint8_t *out, size_t out_cap);
int32_t ds_largest_fragment(const uint8_t *ptr, size_t len, uint8_t *out, size_t out_cap);
int32_t ds_strip_salts(const uint8_t *ptr, size_t len, uint8_t *out, size_t out_cap);
int32_t ds_neutralize_charges(const uint8_t *ptr, size_t len, uint8_t *out, size_t out_cap);
int32_t ds_normalize_smiles(const uint8_t *ptr, size_t len, uint8_t *out, size_t out_cap);
int32_t ds_fragment_parent(const uint8_t *ptr, size_t len, uint8_t *out, size_t out_cap);
int32_t ds_mcs_smarts(const uint8_t *a_ptr, size_t a_len,
                      const uint8_t *b_ptr, size_t b_len,
                      uint8_t *out, size_t out_cap);
int32_t ds_mcs_json(const uint8_t *a_ptr, size_t a_len,
                    const uint8_t *b_ptr, size_t b_len,
                    uint8_t *out, size_t out_cap);
int32_t ds_scaffold_network_json(const uint8_t *ptr, size_t len, uint8_t *out, size_t out_cap);
int32_t ds_num_h_acceptors(const uint8_t *ptr, size_t len);
int32_t ds_num_h_donors(const uint8_t *ptr, size_t len);
int32_t ds_num_rotatable_bonds(const uint8_t *ptr, size_t len);
int32_t ds_ring_count(const uint8_t *ptr, size_t len);
int32_t ds_num_aromatic_rings(const uint8_t *ptr, size_t len);
int32_t ds_num_heteroatoms(const uint8_t *ptr, size_t len);
double ds_fraction_csp3(const uint8_t *ptr, size_t len);

// Wildman-Crippen molar refractivity (MolMR). NaN on invalid SMILES.
double ds_mol_mr(const uint8_t *ptr, size_t len);

// Ring-count descriptors (RDKit rdMolDescriptors). -1 on invalid SMILES.
int32_t ds_num_aliphatic_rings(const uint8_t *ptr, size_t len);
int32_t ds_num_saturated_rings(const uint8_t *ptr, size_t len);
int32_t ds_num_aromatic_heterocycles(const uint8_t *ptr, size_t len);
int32_t ds_num_aromatic_carbocycles(const uint8_t *ptr, size_t len);
int32_t ds_num_saturated_heterocycles(const uint8_t *ptr, size_t len);
int32_t ds_num_saturated_carbocycles(const uint8_t *ptr, size_t len);
int32_t ds_num_aliphatic_heterocycles(const uint8_t *ptr, size_t len);
int32_t ds_num_aliphatic_carbocycles(const uint8_t *ptr, size_t len);
int32_t ds_mol_has_substructure(const uint8_t *smiles_ptr, size_t smiles_len,
                                const uint8_t *smarts_ptr, size_t smarts_len);
int32_t ds_mol_substructure_count(const uint8_t *smiles_ptr, size_t smiles_len,
                                  const uint8_t *smarts_ptr, size_t smarts_len);
int32_t ds_mol_substructure_matches_json(const uint8_t *smiles_ptr, size_t smiles_len,
                                         const uint8_t *smarts_ptr, size_t smarts_len,
                                         uint8_t *out, size_t out_cap);

// SMILES with explicit H atoms (verbose bracket form). Returns length written, or -1.
int32_t ds_add_hydrogens(const uint8_t *ptr, size_t len, uint8_t *out, size_t out_cap);

// Morgan/ECFP fingerprint. Writes ceil(n_bits/8) bytes. Returns bytes written or -1.
int32_t ds_morgan_fp_bits(const uint8_t *ptr, size_t len,
                          uint32_t radius, uint32_t n_bits,
                          uint8_t *out, size_t out_cap);

// MACCS keys fingerprint. Writes a fixed 21 bytes (167-bit vector; bits 1..=166
// are the public keys). Returns bytes written (21) or -1.
int32_t ds_maccs_keys(const uint8_t *ptr, size_t len,
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
int32_t ds_mol_block_property(const uint8_t *data, size_t len,
                              const uint8_t *key, size_t key_len,
                              uint8_t *out, size_t cap);
int32_t ds_mol_block_properties_json(const uint8_t *data, size_t len, uint8_t *out, size_t cap);
int32_t ds_mol_block_atoms_json(const uint8_t *data, size_t len, uint8_t *out, size_t cap);
int32_t ds_mol_block_bonds_json(const uint8_t *data, size_t len, uint8_t *out, size_t cap);
int32_t ds_mol_block_json(const uint8_t *data, size_t len, uint8_t *out, size_t cap);
int32_t ds_sdf_property(const uint8_t *data, size_t len, int32_t record_index,
                        const uint8_t *key, size_t key_len,
                        uint8_t *out, size_t cap);
int32_t ds_sdf_properties_json(const uint8_t *data, size_t len, uint8_t *out, size_t cap);
int32_t ds_mol_block_has_3d(const uint8_t *data, size_t len);
double ds_mol_block_centroid_x(const uint8_t *data, size_t len);
double ds_mol_block_centroid_y(const uint8_t *data, size_t len);
double ds_mol_block_centroid_z(const uint8_t *data, size_t len);
double ds_mol_block_radius_of_gyration(const uint8_t *data, size_t len);
double ds_mol_block_min_x(const uint8_t *data, size_t len);
double ds_mol_block_max_x(const uint8_t *data, size_t len);
double ds_mol_block_min_y(const uint8_t *data, size_t len);
double ds_mol_block_max_y(const uint8_t *data, size_t len);
double ds_mol_block_min_z(const uint8_t *data, size_t len);
double ds_mol_block_max_z(const uint8_t *data, size_t len);

// ===================== PDB crate =====================

// format: 0=auto, 1=pdb, 2=cif, 3=xyz
int32_t ds_structure_atom_count(const uint8_t *data, size_t len, uint8_t format);
int32_t ds_structure_chain_count(const uint8_t *data, size_t len, uint8_t format);
int32_t ds_structure_residue_count(const uint8_t *data, size_t len, uint8_t format);
int32_t ds_structure_model_count(const uint8_t *data, size_t len, uint8_t format);
double ds_structure_centroid_x(const uint8_t *data, size_t len, uint8_t format);
double ds_structure_centroid_y(const uint8_t *data, size_t len, uint8_t format);
double ds_structure_centroid_z(const uint8_t *data, size_t len, uint8_t format);
double ds_structure_radius_of_gyration(const uint8_t *data, size_t len, uint8_t format);
double ds_structure_min_x(const uint8_t *data, size_t len, uint8_t format);
double ds_structure_max_x(const uint8_t *data, size_t len, uint8_t format);
double ds_structure_min_y(const uint8_t *data, size_t len, uint8_t format);
double ds_structure_max_y(const uint8_t *data, size_t len, uint8_t format);
double ds_structure_min_z(const uint8_t *data, size_t len, uint8_t format);
double ds_structure_max_z(const uint8_t *data, size_t len, uint8_t format);

// ===================== SELFIES crate =====================

// SMILES → SELFIES. Returns length written, or -1.
int32_t ds_smiles_to_selfies(const uint8_t *ptr, size_t len, uint8_t *out, size_t cap);

// SELFIES → SMILES. Returns length written, or -1.
int32_t ds_selfies_to_smiles(const uint8_t *ptr, size_t len, uint8_t *out, size_t cap);

// Validate SELFIES: returns 1 if valid, 0 otherwise.
int32_t ds_selfies_is_valid(const uint8_t *ptr, size_t len);

} // extern "C"
