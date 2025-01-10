"""
Classifies: CHEBI:32877 primary amine
"""
from rdkit import Chem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule is a primary amine based on its SMILES string.
    A primary amine is characterized by the presence of the -NH2 group attached to a hydrocarbyl group.

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

    # Check for primary amine pattern: a nitrogen with 2 hydrogens
    # attached to either aliphatic or aromatic carbon
    primary_amine_pattern = Chem.MolFromSmarts("[NX3H2][C,c]")  
    # NX3H2: nitrogen with at least two hydrogens
    # [C,c]: attached to carbon, including both aliphatic and aromatic

    if mol.HasSubstructMatch(primary_amine_pattern):
        return True, "Contains primary amine group (-NH2) attached to hydrocarbyl group"

    return False, "Primary amine group not found"