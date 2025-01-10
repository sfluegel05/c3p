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
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Enhanced SMARTS pattern for detecting a lauric acid moiety within an ester
    # This pattern matches a C(=O)O linked to a 12-carbon chain (lauric acid)
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    lauric_acid_pattern = Chem.MolFromSmarts("[CH3]-[CH2]11[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]")  # 12-carbon chain

    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkage found"

    if not mol.HasSubstructMatch(lauric_acid_pattern):
        return False, "No continuous 12-carbon chain (lauric acid) found"

    return True, "Contains a dodecanoate ester group"