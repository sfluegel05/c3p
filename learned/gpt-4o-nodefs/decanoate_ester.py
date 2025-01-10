"""
Classifies: CHEBI:87658 decanoate ester
"""
from rdkit import Chem

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    A decanoate ester has a 10-carbon alkyl chain connected via an ester linkage.

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
    
    # Check for the presence of a decanoate ester substructure
    decanoate_ester_pattern = Chem.MolFromSmarts("CCCCCCCCCC(=O)O[R]")
    
    if not mol.HasSubstructMatch(decanoate_ester_pattern):
        return False, "No decanoate ester pattern found"
    
    # Count of components matching the pattern (i.e., number of decanoate esters)
    match_count = mol.GetSubstructMatches(decanoate_ester_pattern)
    
    if len(match_count) == 0:
        return False, "Decanoate ester pattern is incomplete or missing"
    
    return True, "Contains valid decanoate ester structure"