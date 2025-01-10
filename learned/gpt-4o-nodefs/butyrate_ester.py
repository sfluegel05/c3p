"""
Classifies: CHEBI:50477 butyrate ester
"""
from rdkit import Chem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule is a butyrate ester based on its SMILES string.
    A butyrate ester contains a butyric acid moiety attached via an ester linkage,
    allowing for reasonable structural isomers or branching.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a butyrate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for ester group (C(=O)O) anywhere in the structure
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester functional group found"

    # Check for flexible butyrate group (4C in series, possibly branched)
    butyrate_group_pattern = Chem.MolFromSmarts("C-C-C-C(=O)O")  
    matches = mol.GetSubstructMatches(butyrate_group_pattern)

    if not matches:
        # Check for branched or equivalent forms
        alternative_butyrate_pattern = Chem.MolFromSmarts("[R]-C(=O)O")
        matches = mol.GetSubstructMatches(alternative_butyrate_pattern)
        if not matches:
            return False, "No butyrate ester group found, even in variations"

    return True, "Contains butyrate ester functional group with valid structure"