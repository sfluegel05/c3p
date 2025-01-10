"""
Classifies: CHEBI:47909 3-oxo-Delta(4) steroid
"""
from rdkit import Chem

def is_3_oxo_Delta_4__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(4) steroid based on its SMILES string. 
    A 3-oxo-Delta(4) steroid has a steroid backbone with a 3-oxo group and a Delta(4) double bond.

    Args:
        smiles (str): SMILES string of the chemical entity.

    Returns:
        bool: True if the molecule is a 3-oxo-Delta(4) steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Adjusted and more flexible substructure pattern for identifying a steroid backbone
    steroid_patterns = [
        Chem.MolFromSmarts('C1CC[C@H]2[C@@H](CCC3=CC(=O)CC[C@]23C)[C@H]1CC=C(C)C=O'),  # Recognizes steroid core with flexible attachments
        Chem.MolFromSmarts('C1CCC2C3C(CC4=CC(=O)CCC34)C2C1'),  # Original pattern for comparison
    ]

    steroid_backbone_found = any(mol.HasSubstructMatch(pattern) for pattern in steroid_patterns)
    if not steroid_backbone_found:
        return False, "No steroid backbone found"
        
    # Identify the presence of a 3-oxo group at the C3 position
    oxo_pattern = Chem.MolFromSmarts('C(=O)C')  # Simplified pattern to detect ketone functional groups
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "No 3-oxo group found at the C3 position"
    
    # Check for the presence of a Delta(4) double bond, specifically between Carbon 4 and 5 within the steroid nucleus
    delta_4_pattern = Chem.MolFromSmarts('C=CC=C')  # Approximate pattern for unsaturation between C4 & C5
    if not mol.HasSubstructMatch(delta_4_pattern):
        return False, "No Delta(4) double bond found"
    
    return True, "Molecule classified as 3-oxo-Delta(4) steroid with appropriate features"