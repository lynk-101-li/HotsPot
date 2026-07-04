"""
N-glycosylation risk assessment for binder–target complexes.

Evaluates whether binder residues sterically clash with known or predicted
N-glycosylation sites on the target protein.

Risk tiers (distance to glycosylation Asn Cβ):
  < 8 Å   → HIGH   — direct collision with GlcNAc core likely
  8–15 Å  → MEDIUM — glycan side chains may reach
  > 15 Å  → SAFE
"""

from __future__ import annotations
import numpy as np
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


HIGH_RISK_CUTOFF = 8.0
MEDIUM_RISK_CUTOFF = 15.0


@dataclass
class GlycanSite:
    resseq: int
    resname: str = "ASN"
    motif: str = ""          # e.g. "NxS/T" sequon
    chain: str = "A"


@dataclass
class GlycanRiskResult:
    site: GlycanSite
    min_dist_A: float
    nearest_binder_res: int
    risk: str                # "HIGH" | "MEDIUM" | "SAFE"


def _classify_risk(dist: float) -> str:
    if dist < HIGH_RISK_CUTOFF:
        return "HIGH"
    if dist < MEDIUM_RISK_CUTOFF:
        return "MEDIUM"
    return "SAFE"


def assess_glycan_risk(
    pdb_path: Path,
    glycan_sites: Sequence[GlycanSite],
    chain_target: str = "A",
    chain_binder: str = "B",
) -> list[GlycanRiskResult]:
    """
    For each known glycosylation site, compute the minimum distance
    from any binder atom to the glycan Asn Cβ (proxy for glycan attachment).
    """
    from interface_analysis import parse_pdb_atoms

    atoms = parse_pdb_atoms(pdb_path)
    binder_atoms = [a for a in atoms if a.chain == chain_binder]

    results = []
    for site in glycan_sites:
        asn_cb = next(
            (a for a in atoms
             if a.chain == chain_target
             and a.resseq == site.resseq
             and a.name == "CB"),
            None,
        )
        if asn_cb is None:
            continue

        min_dist = float("inf")
        nearest_res = -1
        for b in binder_atoms:
            d = float(np.linalg.norm(b.coords - asn_cb.coords))
            if d < min_dist:
                min_dist = d
                nearest_res = b.resseq

        results.append(GlycanRiskResult(
            site=site,
            min_dist_A=round(min_dist, 2),
            nearest_binder_res=nearest_res,
            risk=_classify_risk(min_dist),
        ))

    return results


# Default glycosylation sites for the anonymized immune-checkpoint example.
IMMUNE_CHECKPOINT_GLYCAN_SITES = [
    GlycanSite(resseq=64, motif="Site 1", chain="A"),
    GlycanSite(resseq=77, motif="Site 2", chain="A"),
]
