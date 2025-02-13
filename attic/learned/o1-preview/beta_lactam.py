"""
Classifies: CHEBI:35627 beta-lactam
"""
"""
Classifies: beta-lactam
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

    # Define beta-lactam substructure using SMILES of azetidin-2-one
    beta_lactam_smiles = 'C1CNC1=O'
    beta_lactam_pattern = Chem.MolFromSmiles(beta_lactam_smiles)
    if beta_lactam_pattern is None:
        return None, "Invalid beta-lactam pattern"

    # Search for beta-lactam ring
    if mol.HasSubstructMatch(beta_lactam_pattern):
        return True, "Contains beta-lactam ring (four-membered lactam)"
    else:
        return False, "Does not contain beta-lactam ring"