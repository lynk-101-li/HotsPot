# HotsPot

<p align="center">
  <img src="assets/icon.png" width="200" alt="HotsPot icon"/>
</p>

> **Hot**spot residues on antigen–antibody binding epitopes — served fresh, like a hotpot.

HotsPot is an automated computational pipeline for identifying and ranking **hotspot residues** at antigen–antibody binding interfaces. Given a complex structure (or an ensemble), it scores each epitope residue by its energetic and geometric contribution to binding, and returns a ranked shortlist for experimental follow-up.

---

## Concept

A "hotspot" residue is one whose alanine-scanning free energy penalty (ΔΔG) exceeds ~1 kcal/mol — the handful of interface residues that drive the bulk of binding affinity. Pinpointing them accelerates:

- Epitope mapping & vaccine design
- Antibody engineering (affinity maturation, CDR grafting)
- Escape-mutation prediction

---

## Planned Pipeline

```
Input PDB / AlphaFold model
        │
        ▼
  Interface detection
  (SASA, contact distance)
        │
        ▼
  Per-residue scoring
  (Rosetta ΔΔG / FoldX / ML-based)
        │
        ▼
  Hotspot ranking & filtering
        │
        ▼
  Output: ranked residue table + structure annotation
```

---

## Repository Layout

```
HotsPot/
├── src/
│   ├── main.py               # CLI entry point
│   ├── interface_analysis.py # heavy-atom contact & H-bond detection
│   ├── glycan_risk.py        # N-glycosylation clash assessment
│   └── ranking.py            # AF3 metric-based candidate ranking
├── examples/
│   └── immune_checkpoint_case_study/
│                               # miniprotein binder screening example
├── assets/                   # images and visual assets
├── config/                   # scoring weights, distance cutoffs
└── docs/
```

---

## Case Studies

### Immune checkpoint target — miniprotein binder screening

This case study applies HotsPot to rank miniprotein binder candidates for an immune checkpoint target by interface quality, AF3 confidence, and glycosylation clash risk before proceeding to wet-lab validation (BLI → CD → SPR).

See [`examples/immune_checkpoint_case_study/`](examples/immune_checkpoint_case_study/) for the full analysis report and anonymized metrics.

---

## Status

Active development. Core interface analysis and glycan risk modules are functional; alanine-scanning ΔΔG integration in progress.
