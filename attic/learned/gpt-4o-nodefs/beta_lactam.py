"""
Classifies: CHEBI:35627 beta-lactam
"""
from rdkit import Chem

def is_beta_lactam(smiles: str):
    """
    Determines if a molecule is a beta-lactam based on its SMILES string.
    A beta-lactam is characterized by a four-membered lactam ring.

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

    # Define beta-lactam ring pattern (4-membered ring with a nitrogen)
    beta_lactam_pattern = Chem.MolFromSmarts("C1[NR]C(=O)C1")
    
    if mol.HasSubstructMatch(beta_lactam_pattern):
        return True, "Contains the characteristic beta-lactam ring structure"
    else:
        return False, "Does not contain the beta-lactam ring structure"