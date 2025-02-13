"""
Classifies: CHEBI:46722 carbonate ester
"""
from rdkit import Chem

def is_carbonate_ester(smiles: str):
    """
    Determines if a molecule is a carbonate ester based on its SMILES string.
    A carbonate ester contains the O=C(OX)OX motif, where variations may form cyclic or acyclic esters.

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
    
    # General carbonate ester pattern that includes both acyclic and cyclic structures
    general_pattern = Chem.MolFromSmarts("C(=O)(O*)O*")

    # Check if the molecule matches the pattern
    if mol.HasSubstructMatch(general_pattern):
        return True, "Contains carbonate ester structure"
    
    return False, "No carbonate ester pattern found"

# Example usage
# print(is_carbonate_ester("COC(=O)OC"))  # Example for dimethyl carbonate