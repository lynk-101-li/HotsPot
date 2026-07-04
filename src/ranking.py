"""
Candidate ranking from AF3 confidence metrics and interface analysis.
"""

from __future__ import annotations
import pandas as pd
from pathlib import Path


IPTM_THRESHOLD = 0.80
PLDDT_MIN_THRESHOLD = 70.0
PAE_MAX_THRESHOLD = 2.5


def rank_candidates(metrics_csv: Path) -> pd.DataFrame:
    """
    Load a metrics table and return candidates sorted by composite rank.

    Required columns: candidate, iptm, plddt_mean, plddt_min, pae_min_AtoB,
                      pae_min_BtoA, contact_pairs, hbond_pairs
    """
    df = pd.read_csv(metrics_csv)

    # Hard filters
    df = df[df["iptm"] >= IPTM_THRESHOLD]
    df = df[df["plddt_min"] >= PLDDT_MIN_THRESHOLD]
    df = df[df["pae_min_AtoB"] <= PAE_MAX_THRESHOLD]

    # Composite score (higher = better)
    df = df.copy()
    df["composite_score"] = (
        0.40 * df["iptm"]
        + 0.25 * (df["plddt_mean"] / 100)
        + 0.20 * (1 - df["pae_min_AtoB"] / PAE_MAX_THRESHOLD)
        + 0.15 * (df["hbond_pairs"] / df["hbond_pairs"].max())
    )

    return df.sort_values("composite_score", ascending=False).reset_index(drop=True)
