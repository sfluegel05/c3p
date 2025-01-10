"""
Classifies: CHEBI:35381 monosaccharide
"""
from rdkit import Chem

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a monosaccharide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monosaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # At this complexity, let's assume a simple check based on common aspects of monosaccharides
    # Look for multiple hydroxyl groups
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O' and any(bond.GetBondType() == Chem.rdchem.BondType.SINGLE for bond in atom.GetBonds()))
    if hydroxyl_count < 2:
        return False, f"Insufficient hydroxyl groups, found {hydroxyl_count}"

    # Check for ring structure
    ring_info = mol.GetRingInfo()
    if not ring_info.IsInitialized() or ring_info.NumRings() == 0:
        return False, "No ring structure detected"

    # Check for an aldehyde (C=O with terminal H) or ketone group (C=O in the middle)
    has_carbonyl = any(atom.GetSymbol() == 'C' and any(neighbor.GetSymbol() == 'O' and neighbor.GetImplicitValence() == 1 for neighbor in atom.GetNeighbors()) 
                       for atom in mol.GetAtoms())
    if not has_carbonyl:
        return False, "No carbonyl group detected (aldehyde or ketone)"

    return True, "Contains multiple hydroxyl groups, a ring structure, and a carbonyl group characteristic of monosaccharides"