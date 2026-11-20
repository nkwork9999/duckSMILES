"""Batch descriptor reads on a host-provisioned duckSMILES connection."""
from dataclasses import dataclass

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
