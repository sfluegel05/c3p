"""
Classifies: CHEBI:35627 beta-lactam
"""
from rdkit import Chem

def is_beta_lactam(smiles: str):
    """
    Determines if a molecule is a beta-lactam based on its SMILES string.
    
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

    # Define the SMARTS pattern for a beta-lactam ring
    beta_lactam_pattern = Chem.MolFromSmarts("C1CN(C(=O))C1")
    if mol.HasSubstructMatch(beta_lactam_pattern):
        return True, "Contains a 4-membered beta-lactam ring"
    
    return False, "Does not contain a 4-membered beta-lactam ring"