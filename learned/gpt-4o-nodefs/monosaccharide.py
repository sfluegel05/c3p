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

    # Define SMARTS patterns for furanose and pyranose ring structures
    furanose_pattern = Chem.MolFromSmarts('C1OC(O)C(O)C1')  # 5-membered ring
    pyranose_pattern = Chem.MolFromSmarts('C1OC(O)C(O)C(O)C1')  # 6-membered ring

    # Define for open-chain form (example for hexose)
    open_chain_pattern = Chem.MolFromSmarts('O[C@@H](O)[C@H](O)[C@H](O)C=O')

    # Use SMARTS to check for typical monosaccharide ring structures or open-chain form
    if not (mol.HasSubstructMatch(furanose_pattern) or mol.HasSubstructMatch(pyranose_pattern) or mol.HasSubstructMatch(open_chain_pattern)):
        return False, "No ring structure detected (furanose, pyranose, or open-chain)"

    # Check for carbonyl group variants (aldehyde or ketone)
    carbonyl_variants = [Chem.MolFromSmarts('C(=O)[H]'), Chem.MolFromSmarts('C(=O)[C]')]
    if not any(mol.HasSubstructMatch(pattern) for pattern in carbonyl_variants):
        return False, "No carbonyl group detected (aldehyde or ketone)"

    # Check for stereocenters indicating chirality
    chiral_centers = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    if chiral_centers < 1:
        return False, "No chiral centers detected, unlikely to be a monosaccharide"

    return True, "Contains hydroxyl groups, a ring structure or open-chain form, typical carbonyl features, and chiral centers characteristic of monosaccharides"