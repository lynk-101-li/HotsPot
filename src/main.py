"""
HotsPot — hotspot residue screening for binder–target complexes.

Usage:
    python src/main.py --pdb complex.pdb --metrics metrics.csv --output results/
"""

import argparse
from pathlib import Path

from interface_analysis import analyze_interface
from glycan_risk import assess_glycan_risk, IMMUNE_CHECKPOINT_GLYCAN_SITES
from ranking import rank_candidates


def main():
    parser = argparse.ArgumentParser(description="HotsPot screening pipeline")
    parser.add_argument("--pdb", type=Path, help="Complex PDB/CIF file")
    parser.add_argument("--metrics", type=Path, help="AF3 metrics CSV")
    parser.add_argument("--output", type=Path, default=Path("results"))
    parser.add_argument("--chain-target", default="A")
    parser.add_argument("--chain-binder", default="B")
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)

    if args.pdb:
        print(f"Analyzing interface: {args.pdb}")
        report = analyze_interface(args.pdb, args.chain_target, args.chain_binder)
        print(f"  Contact pairs:     {report.contact_pairs}")
        print(f"  H-bond pairs:      {report.hbond_pairs}")
        print(f"  Hydrophobic pairs: {report.hydrophobic_pairs}")

        print("\nGlycan risk (immune-checkpoint example sites):")
        risks = assess_glycan_risk(
            args.pdb, IMMUNE_CHECKPOINT_GLYCAN_SITES, args.chain_target, args.chain_binder
        )
        for r in risks:
            print(f"  N{r.site.resseq} ({r.site.motif}): {r.min_dist_A:.2f} Å → {r.risk}")

    if args.metrics:
        print(f"\nRanking candidates from: {args.metrics}")
        ranked = rank_candidates(args.metrics)
        out = args.output / "ranked_candidates.csv"
        ranked.to_csv(out, index=False)
        print(ranked[["candidate", "iptm", "plddt_mean", "composite_score"]].to_string())
        print(f"\nSaved to {out}")


if __name__ == "__main__":
    main()
