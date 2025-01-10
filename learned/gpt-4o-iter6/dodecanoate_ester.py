"""
Classifies: CHEBI:87659 dodecanoate ester
"""
from rdkit import Chem

def is_dodecanoate_ester(smiles: str):
    """
    Determines if a molecule is a dodecanoate ester based on its SMILES string.
    A dodecanoate ester is a fatty acid ester in which the carboxylic acid component is lauric acid.

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
    
    # Define more contextualized lauroyl ester pattern using SMARTS
    # Match pattens with exactly 12-carbon chain followed by an ester functional group
    lauroyl_pattern_extended = Chem.MolFromSmarts("C(=O)OCCCCCCCCCCCC")

    # Check for the presence of the lauroyl ester pattern
    if mol.HasSubstructMatch(lauroyl_pattern_extended):
        return True, "Contains characteristic lauroyl ester group indicative of dodecanoate ester"
    
    return False, "Does not contain characteristic lauroyl ester group"