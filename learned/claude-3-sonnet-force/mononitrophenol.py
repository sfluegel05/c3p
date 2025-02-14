"""
Classifies: CHEBI:39362 mononitrophenol
"""
"""
Classifies: CHEBI:16581 mononitrophenol

A mononitrophenol is a phenol carrying a single nitro substituent at an unspecified position.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_mononitrophenol(smiles: str):
    """
    Determines if a molecule is a mononitrophenol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mononitrophenol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get rings and check for aromaticity
    rings = mol.GetRingInfo().AtomRings()
    aromatic_rings = [ring for ring in rings if mol.GetAtomWithIdx(ring[0]).GetIsAromatic()]
    if not aromatic_rings:
        return False, "No aromatic rings found"

    # Check for phenol ring
    phenol_rings = []
    for ring in aromatic_rings:
        atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        if (
            len(ring) == 6
            and sum(1 for atom in atoms if atom.GetAtomicNum() == 6) == 6
            and sum(1 for atom in atoms if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 1) == 1
        ):
            phenol_rings.append(ring)

    if not phenol_rings:
        return False, "No phenol ring found"

    # Check for single nitro group
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    if len(nitro_matches) != 1:
        return False, f"Found {len(nitro_matches)} nitro groups, expected 1"

    # Check if nitro group is attached to phenol ring
    nitro_atom_idx = nitro_matches[0][0]
    nitro_ring_idx = mol.GetAtomWithIdx(nitro_atom_idx).GetIsAromatic()
    if nitro_ring_idx is None or nitro_ring_idx not in [ring[0] for ring in phenol_rings]:
        return False, "Nitro group not attached to phenol ring"

    # Check for other substituents on phenol ring
    for ring in phenol_rings:
        atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        if sum(1 for atom in atoms if atom.GetAtomicNum() not in [6, 8, 7]) > 0:
            return False, "Phenol ring has additional substituents other than nitro group"

    return True, "Contains a phenol ring with a single nitro substituent"