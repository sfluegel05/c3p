"""
Classifies: CHEBI:27933 beta-lactam antibiotic
"""
from rdkit import Chem

def is_beta_lactam_antibiotic(smiles: str):
    """
    Determines if a molecule is a beta-lactam antibiotic based on its SMILES string.
    A beta-lactam antibiotic contains a beta-lactam ring (four-membered cyclic amide), 
    typically in a derivative scaffold like penicillins or cephalosporins.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-lactam antibiotic, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Enhanced SMARTS patterns for beta-lactam structures
    beta_lactam_patterns = [
        Chem.MolFromSmarts("C1C(=O)N(C)C1"),  # Core four-membered beta-lactam ring
        Chem.MolFromSmarts("C1C(=O)N(C)C1C"),  # Simple beta-lactam with alkyl substitution
        Chem.MolFromSmarts("C1C(=O)N2C(C)C2S1"),  # Penicillin-like structure
        Chem.MolFromSmarts("C1C2=C(C)SC(N2C(N1)=O)C(=O)"),  # Cephalosporin-like structure
        Chem.MolFromSmarts("C1C(=O)N2C2C(=O)C1N"),  # Monobactam-like structure
    ]
    
    # Check each pattern for a match
    for pattern in beta_lactam_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains a beta-lactam ring"

    return False, "No beta-lactam ring found"