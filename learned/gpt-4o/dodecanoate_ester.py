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

    # Define SMARTS pattern for dodecanoate (lauric acid) ester: O=C(OCCCCCCCCCCC) - 12 C chain
    dodecanoate_pattern = Chem.MolFromSmarts("O=C(OCCCCCCCCCCC)")

    if not mol.HasSubstructMatch(dodecanoate_pattern):
        return False, "No dodecanoate ester group found"

    return True, "Contains a dodecanoate ester group"