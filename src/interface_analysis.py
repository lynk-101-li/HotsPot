"""
Heavy-atom interface analysis from CIF/PDB structures.

Computes contact pairs, H-bond candidates, and hydrophobic contacts
between two chains in a predicted complex.
"""

from __future__ import annotations
import numpy as np
from pathlib import Path
from dataclasses import dataclass, field
from typing import Sequence


CONTACT_CUTOFF = 5.0       # Å — all heavy-atom contacts
HBOND_CUTOFF = 3.5         # Å — polar contacts (potential H-bonds)
HYDROPHOBIC_CUTOFF = 6.5   # Å — CB–CB distance for hydrophobic pairs

HYDROPHOBIC_AA = {"ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TRP", "PRO", "TYR"}
POLAR_ATOMS = {"N", "O", "S"}


@dataclass
class Atom:
    chain: str
    resname: str
    resseq: int
    name: str
    coords: np.ndarray


@dataclass
class InterfaceReport:
    contact_pairs: int = 0
    hbond_pairs: int = 0
    hydrophobic_pairs: int = 0
    binder_contact_residues: int = 0
    target_contact_residues: int = 0
    details: list = field(default_factory=list)


def parse_pdb_atoms(pdb_path: Path) -> list[Atom]:
    atoms = []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            chain = line[21].strip()
            resname = line[17:20].strip()
            resseq = int(line[22:26])
            name = line[12:16].strip()
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            atoms.append(Atom(chain, resname, resseq, name, np.array([x, y, z])))
    return atoms


def analyze_interface(
    pdb_path: Path,
    chain_target: str = "A",
    chain_binder: str = "B",
) -> InterfaceReport:
    """
    Compute inter-chain contact statistics between target and binder chains.
    """
    atoms = parse_pdb_atoms(pdb_path)
    target = [a for a in atoms if a.chain == chain_target]
    binder = [a for a in atoms if a.chain == chain_binder]

    report = InterfaceReport()
    binder_res_set: set[int] = set()
    target_res_set: set[int] = set()

    for a in binder:
        for b in target:
            dist = float(np.linalg.norm(a.coords - b.coords))

            if dist < CONTACT_CUTOFF:
                report.contact_pairs += 1
                binder_res_set.add(a.resseq)
                target_res_set.add(b.resseq)

                if dist < HBOND_CUTOFF and (
                    a.name[0] in POLAR_ATOMS and b.name[0] in POLAR_ATOMS
                ):
                    report.hbond_pairs += 1

            if a.name == "CB" and b.name == "CB" and dist < HYDROPHOBIC_CUTOFF:
                if a.resname in HYDROPHOBIC_AA and b.resname in HYDROPHOBIC_AA:
                    report.hydrophobic_pairs += 1

    report.binder_contact_residues = len(binder_res_set)
    report.target_contact_residues = len(target_res_set)
    return report
