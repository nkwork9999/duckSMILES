# Davies hydrophilic–lipophilic balance

`chem_davies_hlb(counts MAP(VARCHAR, BIGINT)) → DOUBLE` calculates a
dimensionless HLB from manually assigned functional group counts. Equation (9)
adds 7 to the sum of the signed contributions in Table III. Table II reports
rounded comparison values. [J. T. Davies (1957), *A Quantitative Kinetic Theory
of Emulsion Type, I. Physical Chemistry of the Emulsifying Agent*, pp. 426–438](https://es.firp-ula.org/wp-content/uploads/2019/07/57_Chap_Davies.pdf)

The related minichannel study predicts water-in-oil droplet size from flow and
physical property inputs. The DuckSMILES calculation covers the group-count HLB
described above. [Scientific Reports, doi:10.1038/s41598-026-66575-1](https://doi.org/10.1038/s41598-026-66575-1)

## Design

Evaluate `(7000 + Σ(count × signed_milli)) / 1000` with signed 128-bit products
and accumulation. Convert to double only after the integer accumulation.
The production scalar and `chem_davies_hlb_groups()` share one dictionary.
The table exposes `group_id VARCHAR` and `contribution DECIMAL(6,3)`.

| Group ID | Signed contribution in thousandths |
| --- | ---: |
| sodium_sulfate | 38700 |
| potassium_carboxylate | 21100 |
| sodium_carboxylate | 19100 |
| tertiary_amine | 9400 |
| sorbitan_ester | 6800 |
| free_ester | 2400 |
| carboxylic_acid | 2100 |
| free_hydroxyl | 1900 |
| ether | 1300 |
| sorbitan_hydroxyl | 500 |
| hydrocarbon | -475 |
| ethylene_oxide | 330 |
| propylene_oxide | -150 |

Count the hydrocarbon groups listed in Table III under `hydrocarbon`. When
assigning a derived residue, count that residue once and omit the component
groups for the same portion of the molecule.

Keys are case-sensitive. Missing keys have count zero. A NULL map, empty map,
unknown key (including with count zero), NULL or negative count, or all-zero
counts returns SQL NULL for that row. A calculated zero is a valid score.

## SQL examples

```sql
LOAD ducksmiles;
SELECT chem_davies_hlb(MAP {'free_hydroxyl': 1, 'hydrocarbon': 4}); -- 7.0
SELECT * FROM chem_davies_hlb_groups() ORDER BY group_id;
```

[`scripts/davies_hlb.sql`](../scripts/davies_hlb.sql) defines the standalone
reference macro `davies_hlb_ref`. Run it in a separate connection to compare
the reference result with the production API.

## Test specification

- Verify both column types and all 13 contributions, including their signs.
- Verify n-butanol (1 free hydroxyl, 4 hydrocarbons) at 7.0, cetyl alcohol
  (1 free hydroxyl, 16 hydrocarbons) at 1.3, sodium lauryl sulfate
  (1 sodium sulfate, 12 hydrocarbons) at 40, and sodium oleate
  (1 sodium carboxylate, 17 hydrocarbons) at 18.025 (18 in Table II).
- Verify mixed group types including derived residues, reversed map key order,
  omitted and explicit zero counts, and a valid score of zero.
- Verify mixed valid and invalid rows, each NULL rule, filtered inputs, and
  products and cancellation with BIGINT maximum counts.
- Verify catalog registration and metadata after extension load, compare the
  standalone TEMP macro in a separate connection, and run the materials suite.
