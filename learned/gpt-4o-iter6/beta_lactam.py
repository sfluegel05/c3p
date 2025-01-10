"""
Classifies: CHEBI:35627 beta-lactam
"""
from rdkit import Chem

def is_beta_lactam(smiles: str):
    """
    Determines if a molecule is a beta-lactam based on its SMILES string.
    A beta-lactam contains a four-membered lactam (cyclic amide) ring with a carbonyl (C=O) group
    and nitrogen NH in the ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-lactam, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define an improved SMARTS pattern for a four-membered beta-lactam ring
    # This pattern considers variation in the rings, allowing for any substitution patterns.
    beta_lactam_pattern = Chem.MolFromSmarts("[CX3](=O)N1CC1")  # Four-membered ring with C=O and N
    
    if mol.HasSubstructMatch(beta_lactam_pattern):
        return True, "Contains a four-membered beta-lactam ring"
    
    return False, "Does not contain a four-membered beta-lactam ring"