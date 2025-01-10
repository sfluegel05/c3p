"""
Classifies: CHEBI:27933 beta-lactam antibiotic
"""
from rdkit import Chem

def is_beta_lactam_antibiotic(smiles: str):
    """
    Determines if a molecule is a beta-lactam antibiotic based on its SMILES string.
    A beta-lactam antibiotic contains a beta-lactam ring: a 4-membered ring with one nitrogen atom.

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
    
    # Define the SMARTS pattern for beta-lactam ring
    beta_lactam_pattern = Chem.MolFromSmarts("C1CNC(=O)1")
    
    # Check if the molecule contains the beta-lactam ring
    if mol.HasSubstructMatch(beta_lactam_pattern):
        return True, "Contains a beta-lactam ring"
    else:
        return False, "No beta-lactam ring found"