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

    # Natural product core patterns (expanded with more diverse possible cores)
    core_patterns = [
        Chem.MolFromSmarts("C1CCC(C2CCCC3CNC4C3=C2CCCC4)N1"),  # Example steroid-like core (e.g., tropane alkaloid)
        Chem.MolFromSmarts("O=C1[C@H](OC)OC[C@@H](O)[C@@H]1O"),  # Simple lactone (cyclic ester) 
        Chem.MolFromSmarts("CC1[C@H](O)CCC(=O)C1"),  # Ketone moiety in polyketides
        Chem.MolFromSmarts("c1ccccc1"),  # Basic aromatic ring
    ]
    
    # Synthetic modification patterns (expanded with broader range)
    synthetic_patterns = [
        Chem.MolFromSmarts("[C,c][C](=O)[O,N,c]"),  # wider ester/amide patterns
        Chem.MolFromSmarts("[#6](=O)N"),  # Amide linkage
        Chem.MolFromSmarts("Cl"),  # Presence of chlorine indicating likely synthetic modification
        Chem.MolFromSmarts("[#7]C(=O)C"),  # N-1 carbon amide (peptide-like bond)
        Chem.MolFromSmarts("[#8]C(=O)"),  # General ester linkage
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