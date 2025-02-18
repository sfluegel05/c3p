"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_semisynthetic_derivative(smiles: str):
    """
    Determines if a molecule is likely classified as a semisynthetic derivative.
    It attempts to identify natural product-like cores with synthetic structural modifications.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely semisynthetic, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define potential natural product core patterns
    core_patterns = [
        Chem.MolFromSmarts("[#6]1-[#6]=[#6]-[#6]2-[#6]=[#6]-[#6]-[#6]3-[#6]-1-[#6]-[#6]4-[#6]-[#6](#[#6])-[#6]2-[#6]-[#6]3-[#6]-4"), # steroid-like
        Chem.MolFromSmarts("c1ccccc1"),  # Aryl ring (present in many natural products)
        Chem.MolFromSmarts("[#6]1-[#6]-[#8]-[#6]-[#6]=[#6]-[#8]-1"),  # lactone
    ]
    
    # Define synthetic modification patterns
    synthetic_patterns = [
        Chem.MolFromSmarts("[#6][C](=[O,N])[O,N,c]"),  # general ester/amide
        Chem.MolFromSmarts("[#8]C(=O)"),  # carboxylic ester
        Chem.MolFromSmarts("[#6](=O)C"),  # amide linkage
    ]
    
    # Check for natural core patterns
    core_matches = any(mol.HasSubstructMatch(pattern) for pattern in core_patterns)
    
    # Check for synthetic modifications
    synthetic_matches = any(mol.HasSubstructMatch(pattern) for pattern in synthetic_patterns)
    
    # Evaluate matches and make a decision
    if core_matches and synthetic_matches:
        reason = "Molecule has natural product-like core with synthetic modifications"
        return True, reason
    
    return False, "No characteristic semisynthetic derivative patterns detected"