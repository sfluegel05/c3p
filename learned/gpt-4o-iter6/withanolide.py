"""
Classifies: CHEBI:74716 withanolide
"""
from rdkit import Chem

def is_withanolide(smiles: str):
    """
    Determines if a molecule is a withanolide based on its SMILES string.
    
    A withanolide is characterized as a steroid lactone with specific structural features,
    such as a C28 steroid with a modified side chain forming a lactone ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a withanolide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Pattern for a steroid core - prioritize 6-6-6-5 rings
    steroid_core_pattern = Chem.MolFromSmarts("C1CCC2C3CCC4C=C(C)C(CCC4C3CCC12)C")
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "No valid steroid core structure found"
    
    # Pattern for the presence of a lactone ring, typically -C(=O)O- in a 5-membered ring
    lactone_pattern = Chem.MolFromSmarts("O=C1OC=CC1")
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone ring found"
    
    # Ensure sufficient number of carbons for a C28 steroid
    carbon_count = sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    if carbon_count < 26 or carbon_count > 30:
        return False, f"Molecule has {carbon_count} carbons, expected around 28"
    
    # Check for common functionalities like hydroxyl (often seen in withanolides)
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4H2][OH]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group found, typically present in withanolides"
    
    return True, "The molecule contains key features of a withanolide"