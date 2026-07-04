# Case Study: Hotspot Screening for Immune Checkpoint-Targeting Miniprotein Binders

This case study demonstrates the HotsPot pipeline applied to a set of miniprotein candidates targeting an immune checkpoint protein.

## Background

The target is an immune checkpoint protein with relevance to solid-tumor immunotherapy. We used this system to validate the HotsPot screening workflow before proceeding to wet-lab validation.

## Pipeline Overview

```
AF3 Structural Predictions (rank0, 10 recycling steps)
              │
              ▼
    Interface Quality Analysis
    - Heavy-atom contact pairs  (< 5 Å)
    - Polar / H-bond contacts   (< 3.5 Å)
    - Hydrophobic contacts      (CB–CB < 6.5 Å)
              │
              ▼
    AF3 Confidence Metrics
    - ipTM, pTM, rank_score
    - Per-chain pLDDT (mean & min)
    - Interfacial PAE (A→B, B→A)
              │
              ▼
    Glycosylation Risk Assessment
    - N-glycosylation sites on the target protein: Site 1, Site 2
    - Distance thresholds:
        < 8 Å  → high collision risk with GlcNAc core
        8–15 Å → potential clash with glycan side chains
        > 15 Å → safe
              │
              ▼
    Candidate Ranking & Experimental Design
```

## Key Metrics (anonymized candidates)

| Candidate | ipTM | pLDDT (mean/min) | pae_min A→B | Glycan site 1 dist (Å) | Glycan site 2 dist (Å) | Decision |
|-----------|------|-----------------|-------------|------------------------|------------------------|----------|
| Cand-A    | 0.90 | 95.4 / 81.5     | 1.12        | 13.38                  | 17.31                  | Priority 1 |
| Cand-B    | 0.89 | 95.6 / 77.1     | 1.28        | 9.14                   | 14.11                  | Priority 2 |
| Cand-C    | 0.86 | 93.5 / 75.2     | 1.38        | 7.49                   | 18.40                  | Priority 3 (borderline glycan risk) |
| Cand-D    | 0.90 | 94.6 / 79.3     | 1.10        | 18.20                  | 3.48                   | Downgraded (high glycan clash risk) |
| Cand-E    | 0.71 | 82.9 / 51.9     | 3.68        | -                      | -                      | Excluded (low confidence) |

## Experimental Validation Design

**BLI screening** (all candidates, His-tag protein):
- Rapid binding/non-binding readout
- Rough K_D estimation

**CD spectroscopy** (parallel with BLI):
- Confirm α-helical fold (negative peaks at 208/222 nm)

**SPR kinetics** (BLI-positive candidates only, His-tag removed):
- Precise k_a, k_d, K_D determination

## Notes on His-tag Placement

All candidates use a **C-terminal His₆-tag** (pET29b(+) default, NdeI/XhoI). Rationale:
- C-terminus is distal from the binding interface for Priority 1 and 2 candidates
- Consistent with Baker lab validation standards (Watson et al. *Nature* 2023; *Nat Commun* 2024, 2025)
- TEV cleavage site (ENLYFQS) included for SPR-stage tag removal

## References

- Watson et al. "De novo design of protein structure and function with RFdiffusion." *Nature* **620**, 1089–1100 (2023).
- Cao et al. "Design of protein-binding proteins from target structure alone." *Nature* **605**, 551–560 (2022).
- "De novo design of mini-protein binders broadly neutralizing *C. difficile* toxin B variants." *Nat Commun* (2024).
