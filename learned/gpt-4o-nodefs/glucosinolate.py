"""
Classifies: CHEBI:24279 glucosinolate
"""
from rdkit import Chem

def is_glucosinolate(smiles: str):
    """
    Determines if a molecule is a glucosinolate based on its SMILES string.
    Glucosinolates contain a glucose moiety, a sulfur-linked aglycone, and an N-sulfooxyimino group.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glucosinolate, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for glucosinolate components
    glucose_pattern = Chem.MolFromSmarts("C1[C@H]([C@H]([C@H]([C@@H](C1O)O)O)O)O")
    thioglucose_pattern = Chem.MolFromSmarts("SC1[C@H]([C@H]([C@H]([C@@H](C1O)O)O)O)CO")
    sulfooxyimino_group_pattern = Chem.MolFromSmarts("N=OS(O)(=O)=O")

    # Check for glucose moiety
    if not mol.HasSubstructMatch(glucose_pattern):
        return False, "No glucose moiety found"

    # Check for thioglucose linkage
    if not mol.HasSubstructMatch(thioglucose_pattern):
        return False, "Thioglucose linkage not found"

    # Check for N-sulfooxyimino group
    if not mol.HasSubstructMatch(sulfooxyimino_group_pattern):
        return False, "N-sulfooxyimino group not found"

    # If all key components are present, classify as glucosinolate
    return True, "Contains glucose moiety, thioglucose linkage, and N-sulfooxyimino group consistent with glucosinolates"