"""
Classifies: CHEBI:17334 penicillin
"""
from rdkit import Chem

def is_penicillin(smiles: str):
    """
    Determines if a molecule is a penicillin based on its SMILES string.
    Penicillins are characterized by a beta-lactam ring and a thiazolidine ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a penicillin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for beta-lactam and thiazolidine rings
    beta_lactam_pattern = Chem.MolFromSmarts("C1C(=O)N(C)C1")
    thiazolidine_pattern = Chem.MolFromSmarts("C1C(SC)NC1")
    
    # Check for beta-lactam ring
    if not mol.HasSubstructMatch(beta_lactam_pattern):
        return False, "No beta-lactam ring found"
    
    # Check for thiazolidine ring
    if not mol.HasSubstructMatch(thiazolidine_pattern):
        return False, "No thiazolidine ring found"

    return True, "Contains both beta-lactam and thiazolidine rings indicative of penicillins"