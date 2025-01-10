"""
Classifies: CHEBI:87658 decanoate ester
"""
from rdkit import Chem

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    A decanoate ester results from the esterification of decanoic acid with an alcohol or phenol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a decanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for decanoate ester: 10-carbon chain ending in ester linkage
    decanoate_pattern = Chem.MolFromSmarts("CCCCCCCCCC(=O)OC")  # Linear chain with ester linkage

    # Check if the decanoate pattern matches
    if mol.HasSubstructMatch(decanoate_pattern):
        return True, "Contains decanoate ester structure: recognized as decanoate ester"
    
    return False, "No decanoate ester structure found"