"""
Classifies: CHEBI:32877 primary amine
"""
from rdkit import Chem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule is a primary amine based on its SMILES string.
    A primary amine is characterized by the presence of the -NH2 group attached to a hydrocarbyl chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a primary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for primary amine pattern: a nitrogen with 2 hydrogens and a single carbon bond
    primary_amine_pattern = Chem.MolFromSmarts("[CX4][NH2]")  # CX4 is a carbon to which aliphatic hydrogens are added
    if mol.HasSubstructMatch(primary_amine_pattern):
        return True, "Contains primary amine group (-NH2) attached to hydrocarbyl chain"

    return False, "Primary amine group not found"