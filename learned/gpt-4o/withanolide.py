"""
Classifies: CHEBI:74716 withanolide
"""
from rdkit import Chem

def is_withanolide(smiles: str):
    """
    Determines if a molecule is a withanolide based on its SMILES string.
    Withanolides are C28 steroid lactones with modified side chains forming lactone rings and substituted derivatives.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a withanolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (None, "Invalid SMILES string")

    # Check for the steroid backbone (pattern tuned for flexibility)
    steroid_patterns = [
        Chem.MolFromSmarts("C1C2CC3CC4(C)CCC1C4C=C3C2"),  # A more flexible core steroid pattern
    ]
    
    if not any(mol.HasSubstructMatch(pattern) for pattern in steroid_patterns):
        return (False, "No steroid backbone detected")
    
    # Check for lactone ring (flexible size and configurations)
    lactone_patterns = [
        Chem.MolFromSmarts("C1OC(=O)C=CC1"),  # Gamma-lactone
        Chem.MolFromSmarts("C1OC(=O)CC1"),    # Delta-lactone
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in lactone_patterns):
        return (False, "No lactone group found")
    
    # Consider presence of oxygens which may indicate functional groups
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if oxygen_count < 3:  # Assuming minimum functionalization
        return (False, "Insufficient oxygen-containing functionalities for withanolides")

    return (True, "Contains features consistent with withanolides: steroid backbone, lactone ring, and functional groups")