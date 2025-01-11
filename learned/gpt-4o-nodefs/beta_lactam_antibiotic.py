"""
Classifies: CHEBI:27933 beta-lactam antibiotic
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_lactam_antibiotic(smiles: str):
    """
    Determines if a molecule is a beta-lactam antibiotic based on its SMILES string.
    A beta-lactam antibiotic should have the four-membered beta-lactam ring structure.

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

    # Define the beta-lactam ring SMARTS pattern
    beta_lactam_pattern = Chem.MolFromSmarts("C1C(C(=O)N1)")
    
    # Check for the presence of the beta-lactam ring
    if mol.HasSubstructMatch(beta_lactam_pattern):
        return True, "Contains beta-lactam core structure"

    return False, "Does not contain beta-lactam core structure"