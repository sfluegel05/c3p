"""
Classifies: CHEBI:46722 carbonate ester
"""
from rdkit import Chem

def is_carbonate_ester(smiles: str):
    """
    Determines if a molecule is a carbonate ester based on its SMILES string.
    A carbonate ester typically has the functional group O=C(O-)O- attached to organic groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a carbonate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES to get molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for carbonate ester
    # Ensure we only match the carbonate ester functional group and not other ester-like groups
    patterns = [
        Chem.MolFromSmarts("O=C(O)O"),                 # linear carbonate ester
        Chem.MolFromSmarts("O=C1OC1"),                 # cyclic carbonate like ethylene carbonate
        Chem.MolFromSmarts("[$(C(=O)(O)(OC))]")       # general pattern for the tri-oxo carbonate
    ]

    # Check if the molecule matches any known carbonate ester pattern
    for pattern in patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains carbonate ester functional group"
    
    return False, "No carbonate ester functional group found"