"""
Classifies: CHEBI:46722 carbonate ester
"""
from rdkit import Chem

def is_carbonate_ester(smiles: str):
    """
    Determines if a molecule is a carbonate ester based on its SMILES string.
    A carbonate ester has the functional group O=C(O-)O- attached to organic groups.

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

    # Define the SMARTS pattern for a carbonate ester: R-O-C(=O)-O-R'
    # This pattern ensures that C=O is bonded to two O atoms, each bonded to carbon atoms
    carbonate_ester_pattern = Chem.MolFromSmarts("[$([#6])]OC(=O)O[$([#6])]")

    # Check if the molecule matches the carbonate ester pattern
    if mol.HasSubstructMatch(carbonate_ester_pattern):
        return True, "Contains carbonate ester functional group"
    else:
        return False, "No carbonate ester functional group found"