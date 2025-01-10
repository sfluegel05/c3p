"""
Classifies: CHEBI:35627 beta-lactam
"""
from rdkit import Chem

def is_beta_lactam(smiles: str):
    """
    Determines if a molecule is a beta-lactam based on its SMILES string.
    A beta-lactam contains a four-membered lactam (cyclic amide) ring including
    the carbonyl group (C=O) and nitrogen within the ring.

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

    # Define a SMARTS pattern for a four-membered beta-lactam ring
    # Allow for substitutions but recognize the crucial parts: C=O and N in a 4-membered ring
    beta_lactam_pattern = Chem.MolFromSmarts("C1=O-[NH]-C-C1")  # A loose pattern for a beta-lactam ring
    
    # Check for the presence of a beta-lactam structure
    if mol.HasSubstructMatch(beta_lactam_pattern):
        return True, "Contains a four-membered beta-lactam ring"
    
    # If pattern match fails
    return False, "Does not contain a four-membered beta-lactam ring"