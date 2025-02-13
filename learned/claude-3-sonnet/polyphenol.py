"""
Classifies: CHEBI:26195 polyphenol
"""
"""
Classifies: CHEBI:24851 polyphenol

A polyphenol is defined as a member of the class of phenols that contain 2 or more benzene rings
each of which is substituted by at least one hydroxy group.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyphenol(smiles: str):
    """
    Determines if a molecule is a polyphenol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyphenol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify aromatic rings
    aromatic_rings = Chem.GetSymmSSSR(mol)
    num_aromatic_rings = sum(1 for ring in aromatic_rings if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring))
    if num_aromatic_rings < 2:
        return False, "Less than 2 aromatic rings"

    # Define SMARTS pattern for hydroxyl group
    hydroxy_pattern = Chem.MolFromSmarts("O[H]")

    # Check if each aromatic ring has at least one hydroxyl group
    rings_with_hydroxy = set()
    for ring in aromatic_rings:
        ring_atoms = set(ring)
        has_hydroxy = False
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if mol.HasSubstructMatch(hydroxy_pattern, [atom_idx]):
                has_hydroxy = True
                break
            for neighbor_idx in atom.GetNeighbors():
                if neighbor_idx in ring_atoms:
                    continue
                neighbor = mol.GetAtomWithIdx(neighbor_idx)
                if mol.HasSubstructMatch(hydroxy_pattern, [neighbor_idx]):
                    has_hydroxy = True
                    break
        if has_hydroxy:
            rings_with_hydroxy.add(id(tuple(ring)))

    if len(rings_with_hydroxy) < num_aromatic_rings:
        return False, "Not all aromatic rings have a hydroxyl group"

    return True, "Contains 2 or more aromatic rings, each with at least one hydroxyl group"