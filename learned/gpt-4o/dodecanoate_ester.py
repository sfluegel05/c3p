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

    # Enhanced SMARTS pattern for dodecanoate ester
    # Matching a 12-carbon chain (flexible layout) connected via ester linkage
    # Using a more general pattern for the ester bond while identifying the 12-carbon chain
    dodecanoate_pattern = Chem.MolFromSmarts("C(=O)O[CH2,-][CH2,-][CH2,-][CH2,-][CH2,-][CH2,-][CH2,-][CH2,-][CH2,-][CH3,-]")

    if not mol.HasSubstructMatch(dodecanoate_pattern):
        return False, "No dodecanoate ester group found"

    return True, "Contains a dodecanoate ester group"