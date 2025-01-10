"""
Classifies: CHEBI:74716 withanolide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_withanolide(smiles: str):
    """
    Determines if a molecule is a withanolide based on its SMILES string.
    A withanolide typically has a steroid skeleton with a lactone ring and oxygen functionalities.

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
    
    # Look for core steroid structure - ergostane-like
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C3CCC4C5CCC(C1)C2C5CCC3C4")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No core steroid structure found"

    # Look for lactone ring
    lactone_pattern = Chem.MolFromSmarts("C1C=CC(=O)O1")
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone ring found"

    # Check for oxygen functionalities (hydroxyl, ketone)
    oxy_func_patterns = [
        "[OH]",  # Hydroxyl group
        "[CX3](=O)[#6]"  # Ketone
    ]

    found_oxy_func = any(mol.HasSubstructMatch(Chem.MolFromSmarts(patt)) for patt in oxy_func_patterns)
    if not found_oxy_func:
        return False, "Missing expected oxygen functionalities (hydroxyl, ketone)"

    # If matches all the criteria
    return True, "Contains steroid backbone, lactone ring, and oxygen functionalities"