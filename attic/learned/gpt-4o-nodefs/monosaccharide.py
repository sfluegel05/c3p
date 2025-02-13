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

    # Check for multiple hydroxyl groups
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O' and any(bond.GetBondType() == Chem.rdchem.BondType.SINGLE for bond in atom.GetBonds()))
    if hydroxyl_count < 2:
        return False, f"Insufficient hydroxyl groups, found {hydroxyl_count}"

    # Use SMARTS to check for typical monosaccharide ring structures
    furanose_pattern = Chem.MolFromSmarts('C1OC(O)C(O)C1')
    pyranose_pattern = Chem.MolFromSmarts('C1OC(O)C(O)C(O)C1')
    
    if not mol.HasSubstructMatch(furanose_pattern) and not mol.HasSubstructMatch(pyranose_pattern):
        return False, "No ring structure detected (furanose or pyranose)"

    # Check for carbonyl group (aldehyde C=O at end, ketone C=O internal)
    aldehyde_pattern = Chem.MolFromSmarts('C(=O)[H]')
    ketone_pattern = Chem.MolFromSmarts('C(=O)[C]')
    
    if not (mol.HasSubstructMatch(aldehyde_pattern) or mol.HasSubstructMatch(ketone_pattern)):
        return False, "No carbonyl group detected (aldehyde or ketone)"

    # (Optional) Check for stereocenters indicating chirality
    chiral_centers = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    if chiral_centers < 1:
        return False, "No chiral centers detected, unlikely to be a monosaccharide"

    return True, "Contains hydroxyl groups, a ring structure, a carbonyl group, and chiral centers typical of monosaccharides"