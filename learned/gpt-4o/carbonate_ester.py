"""
Classifies: CHEBI:46722 carbonate ester
"""
from rdkit import Chem

def is_carbonate_ester(smiles: str):
    """
    Determines if a molecule is a carbonate ester based on its SMILES string.
    A carbonate ester contains the O=C(OX)OX motif, where variations may form cyclic esters.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbonate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Carbonate ester pattern including acyclic and cyclic forms
    carbonate_ester_patterns = [
        Chem.MolFromSmarts("O=C(O[*])O[*]"),  # Acyclic carbonate ester
        Chem.MolFromSmarts("O=C1OC(=O)O1"),  # Common cyclic carbonate ester e.g., diethylene carbonate
        Chem.MolFromSmarts("O=C(OC)OC")      # Specific ester type e.g., dimethyl carbonate
    ]
    
    # Additional specific patterns for more variations
    additional_patterns = [
        Chem.MolFromSmarts("O=C(O)OC"),     # Monoester variations
        Chem.MolFromSmarts("O=C1OCCOC1"),   # Cyclic with oxygen in ring
    ]

    for pattern in carbonate_ester_patterns + additional_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains carbonate ester structure"

    return False, "No carbonate ester pattern found"

# Example usage
# print(is_carbonate_ester("COC(=O)OC"))  # Example for dimethyl carbonate