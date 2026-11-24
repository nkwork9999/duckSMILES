# Materials, solvent features and liquid viscosity

Seven SQL functions extend the native, dependency-free Rust/C++ extension.
Outputs are versioned JSON strings. SQL NULL means invalid/unsupported input;
JSON null inside a solvent profile means an undefined descriptor. No new runtime
Python, RDKit, pymatgen or chemicals dependency is introduced.

## Composition and perovskite descriptors

`composition_descriptors_json(composition VARCHAR)` accepts positive integer or
fractional element amounts, repeated symbols and nested parentheses (depth <=16).
Examples: `CsPbI3`, `Cs0.5Rb0.5Pb(I0.8Br0.2)3`, `Ca(OH)2`. Input is ASCII,
<=4096 bytes, with total amount <=1,000,000. Decimal points mean fractional
amounts, never hydrate separators. Charges, isotopes, hydrate separators, variable
symbols and implicit MA/FA aliases are unsupported; expand organic cations to
explicit elemental compositions when computing composition-only features.
Amounts or normalized fractions that underflow to zero are rejected rather than
emitting zero occupancies or nonfinite JSON values.

The result contains alphabetically sorted element amounts and atomic fractions,
number of elements, formula amount and formula weight (g/mol of the entered
formula), mean atomic weight, mean/population-standard-deviation atomic number,
ideal mixing entropy divided by R (`-sum(x ln x)`), and stoichiometry Lp norms
for p=2,3,5,7,10. Fractional features are invariant to rescaling a formula; formula
weight and formula amount are not. Atomic weights reuse the pinned RDKit
2026.03.6 table. These are composition descriptors, not crystal structure features.

`chem_perovskite_descriptors(data DOUBLE[][], parameters DOUBLE[])` takes one or
more `[rA,rB,rX,zA,zB,zX]` rows and `[]`. Site charges may be occupancy-weighted
averages. Radii must be finite and positive in one consistent unit; cation
charges must be positive and anion charge negative. Returns original inputs and:

* Goldschmidt `t=(rA+rX)/(sqrt(2)*(rB+rX))`.
* Octahedral `mu=rB/rX`.
* ABX3 charge sum `zA+zB+3*zX` and charge balance (absolute tolerance 1e-10).

A nonzero charge sum is reported rather than silently corrected. Supply radii
appropriate to the oxidation state and coordination number, and retain their
source/definition in your table. No radius table is guessed from an element
symbol. Mixed-site effective radii must be explicitly defined upstream; organic
cation effective radii are not interchangeable with Shannon atomic ionic radii.
No structural relaxation, CIF processing, band gap, synthesis success or phase
stability classification is performed. The two geometric indicators alone do
not establish that a material is a stable perovskite.

```sql
SELECT composition_descriptors_json('Cs0.5Rb0.5Pb(I0.8Br0.2)3');
-- Radii below illustrate the API; store the actual radius source separately.
SELECT chem_perovskite_descriptors([[1.88,1.19,2.20,1.0,2.0,-1.0]], []);
```

## Molecular and mixture features for solvents/additives

`solvent_descriptors_json(smiles VARCHAR)` returns a single connected molecular
profile, including formal charge and these existing native descriptor IDs:

| ID | Meaning / units |
| --- | --- |
| MolWt | Molecular weight, g/mol |
| MolLogP | Wildman-Crippen octanol/water log partition coefficient |
| MolMR | Wildman-Crippen molar refractivity, conventional RDKit units |
| TPSA | Topological polar surface area, square angstrom; RDKit default N/O scope |
| NumHDonors / NumHAcceptors | Hydrogen-bond donor / acceptor counts |
| NumRotatableBonds | Strict rotatable-bond count |
| HeavyAtomCount | Non-hydrogen atom count |
| FractionCSP3 | Fraction of carbon atoms with sp3 hybridization |
| MaxPartialCharge / MinPartialCharge | Gasteiger partial-charge extrema, elementary-charge units |

The calculation reuses workbench preparation and the existing descriptor catalog;
it does not change other APIs' parsing or descriptor definitions. Disconnected
mixtures and molecules without heavy atoms are rejected. Charged additives are
allowed when supported; supply a separate row for each component. The reference
identifier records descriptor definitions, not a guarantee of full RDKit parser
or sanitization equivalence across arbitrary chemistry. See [existing parity
limits](PROFILE_PARITY.md). These descriptors are potential ML inputs, not
viscosity, dielectric constant, Gutmann donor number or experimental solubility.

`chem_mixture_descriptors` takes `[mole_fraction,feature1,feature2,...]` rows and
`[]`, returning weighted means, population standard deviations and ideal mixing
entropy/R. Every feature must be finite; choose and document any missing-value
policy before calling. Fractions must lie in [0,1] and sum to one within 1e-10;
no silent normalization. Weighted descriptors are feature engineering, not a
physical mixing law for the underlying properties. Feature order and units are
those supplied by the caller. Maximum 10,000 components/rows, 128 columns.

```sql
SELECT solvent_descriptors_json('CS(C)=O'); -- DMSO
SELECT solvent_descriptors_json('CN(C)C=O'); -- DMF
WITH components AS (
  SELECT * FROM (VALUES (0.4, 'CS(C)=O'), (0.6, 'CN(C)C=O')) t(x, smiles)
), features AS (
  SELECT x, solvent_descriptors_json(smiles)::JSON AS j FROM components
)
SELECT chem_mixture_descriptors(
  list([x, (j->'descriptors'->>'MolWt')::DOUBLE,
           (j->'descriptors'->>'TPSA')::DOUBLE]), []) FROM features;
```

## Temperature-dependent viscosity, mixing and fitting

All viscosity inputs/outputs are **dynamic viscosity in Pa s** and temperatures
are **K**. `1 cP = 1 mPa s = 0.001 Pa s`. Convert kinematic viscosity using a
measured density first. Pressure, composition and phase must be fixed across a
fit; a mixed-liquid calculation requires component viscosities at the same
specified temperature and pressure. These conditions are caller metadata, not
inferred from numbers. No concentration-independent or shear-rate-independent
prediction is claimed for non-Newtonian precursor inks.

| SQL function | Data rows | Parameters |
| --- | --- | --- |
| `chem_viscosity_arrhenius` | `[temperature_K]` | `[A,B_K]` |
| `chem_viscosity_fit` | `[temperature_K,viscosity_Pa_s]`, at least 3 | `[]` |
| `chem_viscosity_mix` | `[mole_fraction,viscosity_Pa_s]` | `[]` or binary `[G12]` |

The Andrade/Arrhenius form is `ln(eta/(Pa s))=A+B/T`. Evaluation requires
user-supplied coefficients. Fit uses unweighted OLS in **log viscosity** versus
inverse temperature, not least squares in raw viscosity. It reports reusable
`parameters=[A,B]`, `B*R` in J/mol using R=8.31446261815324, training temperature
range, log residuals, RMSE with denominator n, and residual standard deviation
with denominator n-2. Repeated temperatures are allowed if the design has rank;
identical/numerically indistinguishable temperatures are rejected. No monotonic
constraint is imposed: a negative fitted B is reported, not replaced. Check
residuals and validate on held-out temperatures/compositions; evaluation outside
the reported range is extrapolation. Fitting does not demonstrate predictive
accuracy on new solvents or inks.

Mixing uses `ln(eta_mix)=sum(x_i*ln(eta_i))`. For two components with `[G12]`,
adds `x1*x2*G12` (Grunberg-Nissan). G12 is dimensionless, supplied by the user,
and can depend on temperature/composition. `[]` is a zero-interaction assumption,
not a measured precursor-ink model. Zero-fraction components still require
positive finite viscosity. All methods reject invalid dimensions, parameters,
nonfinite inputs and nonpositive temperatures/viscosities; nonrepresentable
results also return SQL NULL. They retain the 10,000-row/128-column bounds.

```sql
SELECT chem_viscosity_fit(
  [[290.0,0.003],[300.0,0.0025],[310.0,0.002]], []);
SELECT chem_viscosity_arrhenius([[298.15],[320.0]], [-12.0,1500.0]);
SELECT chem_viscosity_mix([[0.5,0.001],[0.5,0.004]], []); -- 0.002 Pa s
SELECT chem_viscosity_mix([[0.5,0.001],[0.5,0.004]], [0.8]);
```

## References and validation

* [matminer composition descriptors](https://hackingmaterials.lbl.gov/matminer/matminer.featurizers.composition.html): element fractions and normalized stoichiometry features; no claim of full Magpie feature parity.
* [Perovskite geometric descriptors and their limitations](https://pmc.ncbi.nlm.nih.gov/articles/PMC6368436/).
* [RDKit descriptor definitions](https://www.rdkit.org/docs/GettingStartedInPython.html#descriptor-calculation).
* [Binary Grunberg-Nissan equation and experimental use](https://www.mdpi.com/1420-3049/26/18/5513).
* [Liquid viscosity correlations and mixing methods](https://chemicals.readthedocs.io/chemicals.viscosity.html).

Rust regressions test formulas, invalid inputs, FFI buffers, mixture identities,
fit recovery and numeric failure paths. `test/sql/materials.test` covers the real
SQL boundary. `scripts/check_materials.py` independently compares the loaded
extension against RDKit 2026.03.6, pymatgen composition parsing, scipy regression,
and chemicals logarithmic mixing (development-only Python dependencies).

### Verified checkpoint

The [validation receipt](tasks/evidence/materials-validation.json) records:
1,429 Rust tests passed, zero failed, two pre-existing ignored tests; 4,386 native
SQL/reference assertions passed on macOS arm64 Release, DuckDB 1.5.5. This
includes 15 solvent/additive structures, nine compositions, 100 random radius and
mixture cases, the noisy-temperature regression and selected multi-batch SQL
vectors. Results are numerical/implementation validation, not a measured
viscosity-prediction benchmark. Windows, Linux and Wasm results must be checked
in submission CI separately. Analytics strict Clippy passed; the existing smiles
crate Clippy raw-pointer ABI findings remain outside this change.

Reproduce with a matching built extension:

```sh
CARGO_TARGET_DIR=/private/tmp/ducksmiles-materials-target cargo test --workspace --locked
python scripts/check_materials.py --extension /path/to/ducksmiles.duckdb_extension
```

## Davies hydrophilic–lipophilic balance

```sql
LOAD ducksmiles;
SELECT chem_davies_hlb(MAP {'free_hydroxyl': 1, 'hydrocarbon': 4}); -- 7.0
SELECT * FROM chem_davies_hlb_groups() ORDER BY group_id;
```

Supply manually assigned group counts. The calculation and group identifiers
follow Davies (1957), equation (9) and Table III. [Group assignments, source
and reference SQL](docs/emulsion.md)
