"""
Classifies: CHEBI:35627 beta-lactam
"""
"""
Classifies: Beta-lactam
"""
from rdkit import Chem

def is_beta_lactam(smiles: str):
    """
    Determines if a molecule is a beta-lactam based on its SMILES string.
    A beta-lactam is a lactam in which the amide bond is contained within a four-membered ring,
    which includes the amide nitrogen and the carbonyl carbon.

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

    # Define beta-lactam SMARTS pattern
    beta_lactam_smarts = '[N;R][C;R](=O)[C;R][C;R]'
    beta_lactam_pattern = Chem.MolFromSmarts(beta_lactam_smarts)

    # Check for beta-lactam ring
    if mol.HasSubstructMatch(beta_lactam_pattern):
        return True, "Molecule contains a beta-lactam ring"
    else:
        return False, "No beta-lactam ring found"