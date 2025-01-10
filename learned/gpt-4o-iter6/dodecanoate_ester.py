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

    # Define the full lauroyl ester pattern using SMARTS
    # - Look for a 12-carbon chain that ends with a carbonyl oxygen linkage (-COO)
    lauroyl_ester_pattern = Chem.MolFromSmarts("CCCCCCCCCCCC(=O)O")

    # Check for the presence of the lauroyl ester pattern
    if mol.HasSubstructMatch(lauroyl_ester_pattern):
        return True, "Contains characteristic lauroyl ester group indicative of dodecanoate ester"
    
    return False, "Does not contain characteristic lauroyl ester group"