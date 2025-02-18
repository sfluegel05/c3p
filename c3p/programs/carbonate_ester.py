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

    # Define the SMARTS pattern for a carbonate species
    # Capture linear and cyclic variations of carbonate groups
    carbonate_ester_patterns = [
        Chem.MolFromSmarts("[$([#6])-]-O-C(=O)-O-[$([#6])]"),  # R-O-C(=O)-O-R', linear
        Chem.MolFromSmarts("O=C1OC[*]C1"),  # Cyclic carbonate, e.g., ethylene carbonate
        Chem.MolFromSmarts("C(=O)(O[*])O[*]")  # General template, possible variations within rings
    ]

    # Check if the molecule matches any known carbonate ester pattern
    for pattern in carbonate_ester_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains carbonate ester functional group"
    
    return False, "No carbonate ester functional group found"