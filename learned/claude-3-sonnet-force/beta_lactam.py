"""
Classifies: CHEBI:35627 beta-lactam
"""
"""
Classifies: CHEBI:35458 beta-lactams
A lactam in which the amide bond is contained within a four-membered ring,
which includes the amide nitrogen and the carbonyl carbon.
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

    # Look for 4-membered ring with N and C=O
    beta_lactam_pattern = Chem.MolFromSmarts("[NR2]1[CR2][CR2][CR2]1=O")
    if mol.HasSubstructMatch(beta_lactam_pattern):
        return True, "Contains a 4-membered ring with an amide nitrogen and a carbonyl carbon (key structural feature of beta-lactams)"
    else:
        return False, "Does not contain the key structural feature of beta-lactams (4-membered ring with N and C=O)"