"""
Classifies: CHEBI:27933 beta-lactam antibiotic
"""
from rdkit import Chem

def is_beta_lactam_antibiotic(smiles: str):
    """
    Determines if a molecule is a beta-lactam antibiotic based on its SMILES string.
    A beta-lactam antibiotic contains a beta-lactam ring (four-membered cyclic amide), 
    often as part of a larger bicyclic structure.

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
    
    # Comprehensive SMARTS pattern for a beta-lactam motif including fused rings
    beta_lactam_patterns = [
        Chem.MolFromSmarts("C1C(=O)NC1"),  # Simple beta-lactam
        Chem.MolFromSmarts("C1C(=O)N2C=CN2C1=O"),  # Penicillins (penam)
        Chem.MolFromSmarts("C1C(=O)N2C=CC2S1=O"),  # Cephalosporins (cephem)
        Chem.MolFromSmarts("C1C2C(=O)NC1C(=O)S2"),  # Carbapenems (penem)
        Chem.MolFromSmarts("N1C(=O)C2CC2C1=O"),  # Monobactams
    ]
    
    # Check each distinct pattern for a match
    for pattern in beta_lactam_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains a beta-lactam ring"

    return False, "No beta-lactam ring found"