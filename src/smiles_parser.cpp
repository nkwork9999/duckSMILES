#include "smiles_parser.hpp"

#include <algorithm>
#include <cctype>
#include <set>
#include <stack>

namespace smiles {

static const std::set<std::string> ORGANIC_SUBSET = {"B", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I"};
static const std::set<char> AROMATIC_ATOMS = {'b', 'c', 'n', 'o', 'p', 's'};

static int default_valence(const std::string &sym, bool aromatic) {
	if (sym == "B") return aromatic ? 2 : 3;
	if (sym == "C") return aromatic ? 3 : 4;
	if (sym == "N") return aromatic ? 2 : 3;
	if (sym == "O") return 2;
	if (sym == "P") return 3;
	if (sym == "S") return 2;
	if (sym == "F" || sym == "Cl" || sym == "Br" || sym == "I") return 1;
	return 0;
}

// Parse bracket atom: [NH3+], [Fe+2], [13C@@H], etc.
// pos points after '[', advances past ']'
static bool parse_bracket_atom(const std::string &smi, size_t &pos, Atom &atom) {
	atom.in_bracket = true;

	// Skip isotope
	while (pos < smi.size() && std::isdigit(smi[pos])) pos++;
	if (pos >= smi.size()) return false;

	// Element symbol
	if (smi[pos] == '*') {
		atom.symbol = "*";
		pos++;
	} else if (std::islower(smi[pos])) {
		atom.aromatic = true;
		std::string sym(1, static_cast<char>(std::toupper(smi[pos])));
		pos++;
		if (pos < smi.size() && std::islower(smi[pos])) { sym += smi[pos]; pos++; }
		atom.symbol = sym;
	} else if (std::isupper(smi[pos])) {
		atom.symbol = std::string(1, smi[pos]);
		pos++;
		if (pos < smi.size() && std::islower(smi[pos])) { atom.symbol += smi[pos]; pos++; }
	} else {
		return false;
	}

	// Skip chirality
	while (pos < smi.size() && smi[pos] == '@') pos++;

	// Explicit H
	if (pos < smi.size() && smi[pos] == 'H') {
		pos++;
		if (pos < smi.size() && std::isdigit(smi[pos])) {
			atom.hydrogen = smi[pos] - '0'; pos++;
		} else {
			atom.hydrogen = 1;
		}
	}

	// Charge
	if (pos < smi.size() && (smi[pos] == '+' || smi[pos] == '-')) {
		int sign = (smi[pos] == '+') ? 1 : -1;
		pos++;
		if (pos < smi.size() && std::isdigit(smi[pos])) {
			atom.charge = sign * (smi[pos] - '0'); pos++;
		} else {
			atom.charge = sign;
			while (pos < smi.size() && smi[pos] == (sign > 0 ? '+' : '-')) {
				atom.charge += sign; pos++;
			}
		}
	}

	if (pos >= smi.size() || smi[pos] != ']') return false;
	pos++;
	return true;
}

Molecule parse(const std::string &smi) {
	Molecule mol;
	if (smi.empty()) return mol;

	size_t pos = 0;
	std::stack<int> branch_stack;
	std::map<int, int> ring_openings;
	std::vector<int> degree; // per-atom bond degree (counts double/triple bonds by order)
	int prev_atom = -1;
	int next_bond_order = 1; // 1=single, 2=double, 3=triple

	while (pos < smi.size()) {
		char c = smi[pos];

		// Branch
		if (c == '(') {
			if (prev_atom < 0) return mol;
			branch_stack.push(prev_atom);
			pos++; continue;
		}
		if (c == ')') {
			if (branch_stack.empty()) return mol;
			prev_atom = branch_stack.top();
			branch_stack.pop();
			next_bond_order = 1;
			pos++; continue;
		}

		// Bond symbols — record bond order for next bond
		if (c == '=' || c == '#' || c == '-' || c == ':' || c == '/' || c == '\\') {
			if (c == '=') next_bond_order = 2;
			else if (c == '#') next_bond_order = 3;
			else next_bond_order = 1;
			pos++; continue;
		}

		// Dot (disconnected fragments)
		if (c == '.') {
			prev_atom = -1;
			pos++; continue;
		}

		// Ring closure
		if (std::isdigit(c) || c == '%') {
			int ring_num;
			if (c == '%') {
				pos++;
				if (pos + 1 >= smi.size() || !std::isdigit(smi[pos]) || !std::isdigit(smi[pos+1])) return mol;
				ring_num = (smi[pos] - '0') * 10 + (smi[pos+1] - '0');
				pos += 2;
			} else {
				ring_num = c - '0';
				pos++;
			}
			if (ring_openings.count(ring_num)) {
				int other = ring_openings[ring_num];
				mol.bond_count++;
				degree[other] += next_bond_order;
				degree[prev_atom] += next_bond_order;
				next_bond_order = 1;
				ring_openings.erase(ring_num);
			} else {
				ring_openings[ring_num] = prev_atom;
			}
			continue;
		}

		// Bracket atom
		if (c == '[') {
			pos++;
			Atom atom;
			if (!parse_bracket_atom(smi, pos, atom)) return mol;
			int idx = static_cast<int>(mol.atoms.size());
			mol.atoms.push_back(atom);
			degree.push_back(0);
			if (prev_atom >= 0) {
				mol.bond_count++;
				degree[prev_atom] += next_bond_order;
				degree[idx] += next_bond_order;
				next_bond_order = 1;
			}
			prev_atom = idx;
			continue;
		}

		// Aromatic atom
		if (AROMATIC_ATOMS.count(c)) {
			Atom atom;
			atom.aromatic = true;
			atom.symbol = std::string(1, static_cast<char>(std::toupper(c)));
			pos++;
			int idx = static_cast<int>(mol.atoms.size());
			mol.atoms.push_back(atom);
			degree.push_back(0);
			if (prev_atom >= 0) {
				mol.bond_count++;
				degree[prev_atom] += next_bond_order;
				degree[idx] += next_bond_order;
				next_bond_order = 1;
			}
			prev_atom = idx;
			continue;
		}

		// Organic subset atom
		if (std::isupper(c)) {
			std::string sym(1, c);
			if (pos + 1 < smi.size() && std::islower(smi[pos+1])) {
				std::string two = sym + smi[pos+1];
				if (ORGANIC_SUBSET.count(two)) {
					sym = two; pos++;
				}
			}
			pos++;
			if (!ORGANIC_SUBSET.count(sym)) return mol;

			Atom atom;
			atom.symbol = sym;
			int idx = static_cast<int>(mol.atoms.size());
			mol.atoms.push_back(atom);
			degree.push_back(0);
			if (prev_atom >= 0) {
				mol.bond_count++;
				degree[prev_atom] += next_bond_order;
				degree[idx] += next_bond_order;
				next_bond_order = 1;
			}
			prev_atom = idx;
			continue;
		}

		// Unknown character → invalid
		return mol;
	}

	if (!branch_stack.empty() || !ring_openings.empty()) return mol;
	if (mol.atoms.empty()) return mol;

	// Compute implicit H for non-bracket organic subset atoms
	for (size_t i = 0; i < mol.atoms.size(); i++) {
		auto &atom = mol.atoms[i];
		if (atom.in_bracket) continue; // bracket atoms have explicit H already

		int val = default_valence(atom.symbol, atom.aromatic);
		int implicit_h = val - degree[i];
		if (implicit_h < 0) implicit_h = 0;
		atom.hydrogen = implicit_h;
	}

	mol.valid = true;
	return mol;
}

// ============================================================================
// Molecule methods
// ============================================================================

std::map<std::string, int> Molecule::element_counts() const {
	std::map<std::string, int> counts;
	for (const auto &atom : atoms) {
		counts[atom.symbol]++;
		if (atom.hydrogen > 0) {
			counts["H"] += atom.hydrogen;
		}
	}
	return counts;
}

std::string Molecule::formula() const {
	auto counts = element_counts();
	if (counts.empty()) return "";

	std::string result;
	auto append = [&](const std::string &elem, int count) {
		result += elem;
		if (count > 1) result += std::to_string(count);
	};

	// Hill system: C first, H second, rest alphabetical
	if (counts.count("C")) {
		append("C", counts["C"]);
		counts.erase("C");
		if (counts.count("H")) {
			append("H", counts["H"]);
			counts.erase("H");
		}
	}

	std::vector<std::pair<std::string, int>> rest(counts.begin(), counts.end());
	std::sort(rest.begin(), rest.end());
	for (size_t i = 0; i < rest.size(); i++) {
		append(rest[i].first, rest[i].second);
	}
	return result;
}

int Molecule::heavy_atom_count() const {
	int count = 0;
	for (const auto &atom : atoms) {
		if (atom.symbol != "H") count++;
	}
	return count;
}

int Molecule::num_bonds() const {
	return bond_count;
}

} // namespace smiles
