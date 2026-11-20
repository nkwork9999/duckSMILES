"""Batch descriptor reads on a host-provisioned duckSMILES connection."""
from dataclasses import dataclass
from sys import float_info

from duckdb import DuckDBPyConnection


@dataclass(frozen=True)
class DescriptorRow:
    formula: str | None
    weight: float | None


class MoleculeReader:
    def __init__(self, connection: DuckDBPyConnection):
        self._db = connection

    def describe(self, molecules: list[str | None]) -> list[DescriptorRow]:
        """Compute each descriptor per distinct SMILES, then restore input order."""
        if not isinstance(molecules, list) or any(
            value is not None and type(value) is not str for value in molecules
        ):
            raise TypeError("molecules must be a list of SMILES strings or None")
        rows = self._db.execute("""
            WITH inputs AS (
                SELECT smiles, ord FROM unnest(?::VARCHAR[]) WITH ORDINALITY AS t(smiles, ord)
            ), properties AS MATERIALIZED (
                SELECT smiles, mol_formula(smiles) AS formula, mol_weight(smiles) AS weight
                FROM (SELECT DISTINCT smiles FROM inputs WHERE smiles IS NOT NULL)
            )
            SELECT formula, weight FROM inputs LEFT JOIN properties USING (smiles) ORDER BY ord
        """, [molecules]).fetchall()
        return [DescriptorRow(*row) for row in rows]


@dataclass(frozen=True)
class DaltonMass:
    value: float

    def __post_init__(self):
        if type(self.value) not in (float, int) or not 0 < self.value <= float_info.max:
            raise ValueError("mass must be a finite positive number of daltons")


@dataclass(frozen=True)
class MassMatch:
    reference_index: int
    peak_index: int
    error_ppm: float


class MassWindowSearch:
    """Trusted embedded operation for neutral masses expressed in daltons."""

    def __init__(self, connection: DuckDBPyConnection):
        self._db = connection

    def match(self, references: tuple[DaltonMass, ...], peaks: tuple[DaltonMass, ...],
              tolerance_ppm: float) -> list[MassMatch]:
        """Return inclusive matches in input order, with zero-based indices."""
        if type(tolerance_ppm) not in (int, float) or not 0 <= tolerance_ppm <= 1_000_000:
            raise ValueError("tolerance must be between 0 and 1,000,000 ppm")
        for masses in (references, peaks):
            if type(masses) is not tuple or any(type(mass) is not DaltonMass for mass in masses):
                raise TypeError("expected a tuple of DaltonMass values")
        rows = self._db.execute("""
            WITH reference_masses AS (
                SELECT mass, ord, mass * (?::DOUBLE / 1000000) AS delta
                FROM unnest(?::DOUBLE[]) WITH ORDINALITY AS t(mass, ord)
            ), windows AS (
                SELECT mass, ord, mass - delta AS lo, mass + delta AS hi FROM reference_masses
            ), peaks AS (
                SELECT mass, ord FROM unnest(?::DOUBLE[]) WITH ORDINALITY AS t(mass, ord)
            )
            SELECT w.ord - 1, p.ord - 1, (p.mass - w.mass) / w.mass * 1000000
            FROM windows w JOIN peaks p ON p.mass BETWEEN w.lo AND w.hi
            ORDER BY w.ord, p.ord
        """, [tolerance_ppm, [m.value for m in references], [m.value for m in peaks]]).fetchall()
        return [MassMatch(*row) for row in rows]
