"""
Classifies: CHEBI:74716 withanolide
"""
from rdkit import Chem

def is_withanolide(smiles: str):
    """
    Determines if a molecule is a withanolide based on its SMILES string.
    A withanolide typically has a steroid-like skeleton with some oxygen-containing
    functionalities and a lactone group or similar structural connection.

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
    
    # A more flexible pattern for steroid core structure
    steroid_pattern = Chem.MolFromSmarts("C1CC2CC3C4CCC(C4)CC3CCC2C1")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No flexible core steroid-like structure found"

    # Look for a lactone or similar ester group pattern
    lactone_patterns = [
        Chem.MolFromSmarts("O=C1OC[C@@H]1"),  # Simple lactone ring
        Chem.MolFromSmarts("C1OC(=O)CC1")     # More inclusive lactone/ester pattern
    ]
    
    found_lactone = any(mol.HasSubstructMatch(patt) for patt in lactone_patterns)
    if not found_lactone:
        return False, "No lactone or ester-type cyclic structure found"
        
    # Check for additional oxygen functionalities
    oxy_func_patterns = [
        Chem.MolFromSmarts("[OH]"),           # Hydroxyl group
        Chem.MolFromSmarts("[CX3](=O)[#6]")   # Ketone
    ]
    
    found_oxy_func = any(mol.HasSubstructMatch(patt) for patt in oxy_func_patterns)
    if not found_oxy_func:
        return False, "Missing expected oxygen functionalities (hydroxyl, ketone)"

    # If matches all the criteria
    return True, "Contains steroid-like structure, lactone/ester ring, and oxygen functionalities"