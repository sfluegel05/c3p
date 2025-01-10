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
    
    # Potential steroid core - typically (6-6-6-5) rings but with more flexibility and rearrangements
    steroid_core_patterns = [
        # A more flexible pattern that accounts for different steroid core arrangements
        Chem.MolFromSmarts("C1CCC2C3CCC4[C,C](C)C[C,c]4C3CCC2C1"),
        Chem.MolFromSmarts("C1=CC2CCC3[C,C](C)C[C,c]3C2CCC1"), # other plausible ring arrangements
    ]
    
    # Check if any steroid core pattern is present
    if not any(mol.HasSubstructMatch(pattern) for pattern in steroid_core_patterns):
        return False, "No valid steroid core structure found"
    
    # Ensure presence of lactone ring
    lactone_patterns = [
        Chem.MolFromSmarts("O=C1OC=C[C,C]1"),   # 5-membered lactone ring
        Chem.MolFromSmarts("O=C2C=C1C=CC=C1C=C2"), # other plausible lactone formations        
    ]
    
    # Check for any lactone pattern
    if not any(mol.HasSubstructMatch(pattern) for pattern in lactone_patterns):
        return False, "No lactone ring found"
    
    # Ensure sufficient number of carbons, typical of a C28 steroid
    carbon_count = sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    if carbon_count < 25 or carbon_count > 32:
        return False, f"Molecule has {carbon_count} carbons, expected around 28"
    
    # Known functionalities like hydroxyl groups (OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4H1][OH]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group found, typically often present in withanolides"
    
    return True, "The molecule contains key features of a withanolide"