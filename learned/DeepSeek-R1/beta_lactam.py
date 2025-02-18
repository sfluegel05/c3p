"""
Classifies: CHEBI:35627 beta-lactam
"""
"""
Classifies: Beta-lactam (CHEBI:35627)
"""
from rdkit import Chem

def is_beta_lactam(smiles: str):
    """
    Determines if a molecule is a beta-lactam based on its SMILES string.
    A beta-lactam contains a four-membered ring with an amide bond (azetidin-2-one).

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

    # Define beta-lactam SMARTS pattern: four-membered ring with amide bond
    beta_lactam_pattern = Chem.MolFromSmarts("[NH]1C(=O)CC1")  # Azetidin-2-one core

    # Check for presence of the beta-lactam ring
    has_lactam = mol.HasSubstructMatch(beta_lactam_pattern)

    if has_lactam:
        return True, "Contains a four-membered lactam ring (beta-lactam)"
    else:
        return False, "No four-membered lactam ring detected"