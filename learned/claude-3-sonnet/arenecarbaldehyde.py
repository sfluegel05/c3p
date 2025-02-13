"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
"""
Classifies: CHEBI:32879 arenecarbaldehyde
Any aldehyde in which the carbonyl group is attached to an aromatic moiety.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an arenecarbaldehyde, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for aldehyde group
    aldehyde_pattern = Chem.MolFromSmarts("[CH2]=O")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    if not aldehyde_matches:
        return False, "No aldehyde group found"

    # Look for aromatic ring attached to the aldehyde carbon
    aromatic_pattern = Chem.MolFromSmarts("a1aa[CH2]=O")
    aromatic_matches = mol.GetSubstructMatches(aromatic_pattern)
    if not aromatic_matches:
        return False, "Aldehyde group not attached to an aromatic ring"

    # Check if the aromatic ring is part of the same molecule
    aromatic_ring_idx = aromatic_matches[0][1]
    aromatic_ring = mol.GetAtomWithIdx(aromatic_ring_idx).GetIsAromatic()
    if not aromatic_ring:
        return False, "Aldehyde attached to a non-aromatic ring"

    return True, "Contains an aldehyde group attached to an aromatic ring"