#pragma once

#include <string>
#include <map>
#include <vector>
#include <cstdint>

namespace smiles {

// Parsed atom from SMILES
struct Atom {
	std::string symbol;   // "C", "N", "O", "Cl", "Br", etc.
	int hydrogen = 0;     // explicit H count (from bracket atoms)
	int charge = 0;
	bool aromatic = false;
	bool in_bracket = false;
};

// Parse result
struct Molecule {
	std::vector<Atom> atoms;
	int bond_count = 0;
	bool valid = false;

	// Count atoms by element symbol (includes implicit H)
	std::map<std::string, int> element_counts() const;

	// Molecular formula string (Hill system: C first, H second, rest alphabetical)
	std::string formula() const;

	// Heavy atom count (non-hydrogen)
	int heavy_atom_count() const;

	// Total bond count
	int num_bonds() const;
};

// Parse a SMILES string into a Molecule. Sets valid=false on failure.
Molecule parse(const std::string &smiles);

} // namespace smiles
