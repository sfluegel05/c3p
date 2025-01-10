"""
Classifies: CHEBI:32878 alkene
"""
from rdkit import Chem

def is_alkene(smiles: str):
    """
    Determines if a molecule is an alkene based on its SMILES string.
    An alkene contains at least one carbon-carbon double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkene, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find carbon-carbon double bond pattern (C=C)
    alkene_pattern = Chem.MolFromSmarts("C=C")
    if mol.HasSubstructMatch(alkene_pattern):
        return True, "Contains at least one C=C double bond"
    else:
        return False, "No C=C double bond found"