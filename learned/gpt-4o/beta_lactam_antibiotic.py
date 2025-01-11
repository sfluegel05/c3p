"""
Classifies: CHEBI:27933 beta-lactam antibiotic
"""
from rdkit import Chem

def is_beta_lactam_antibiotic(smiles: str):
    """
    Determines if a molecule is a beta-lactam antibiotic based on its SMILES string.
    A beta-lactam antibiotic contains a beta-lactam ring, which is a four-membered cyclic amide.

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

    # Define beta-lactam ring pattern (four-membered cyclic amide)
    beta_lactam_pattern = Chem.MolFromSmarts('C1CNC1=O')
    if mol.HasSubstructMatch(beta_lactam_pattern):
        return True, "Beta-lactam ring detected"
    
    return False, "Beta-lactam ring not found"