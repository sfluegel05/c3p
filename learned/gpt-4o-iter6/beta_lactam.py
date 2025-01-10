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

    # Define an improved SMARTS pattern for a four-membered beta-lactam ring:
    # It identifies the specific 4-membered ring with C=O and N.
    beta_lactam_pattern = Chem.MolFromSmarts("C1(=O)NCC1")  # Basic beta-lactam ring structure

    # Allow for any substitutions on the carbon atoms but must be four-membered
    if mol.HasSubstructMatch(beta_lactam_pattern):
        return True, "Contains a four-membered beta-lactam ring"
    
    return False, "Does not contain a four-membered beta-lactam ring"