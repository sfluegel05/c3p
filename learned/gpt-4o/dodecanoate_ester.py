"""
Classifies: CHEBI:87659 dodecanoate ester
"""
from rdkit import Chem

def is_dodecanoate_ester(smiles: str):
    """
    Determines if a molecule is a dodecanoate ester based on its SMILES string.
    A dodecanoate ester contains a lauric acid moiety (12-carbon chain) as the acid component.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dodecanoate ester, False otherwise
        str: Reason for classification or error message
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for detecting an ester linkage: C(=O)O
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkage found"

    # SMARTS pattern for detecting a 12-carbon chain
    lauric_acid_pattern = Chem.MolFromSmarts("C(=O)OCCCCCCCCCCCC")
    if not mol.HasSubstructMatch(lauric_acid_pattern):
        return False, "No continuous 12-carbon chain (lauric acid) found as part of the ester"
    
    return True, "Contains a dodecanoate ester group"