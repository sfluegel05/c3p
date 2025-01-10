"""
Classifies: CHEBI:24279 glucosinolate
"""
from rdkit import Chem

def is_glucosinolate(smiles: str):
    """
    Determines if a molecule is a glucosinolate based on its SMILES string.
    Glucosinolates contain either a glucose-derived moiety, a sulfur linkage,
    and an N-sulfooxyimino group.
    
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

    # Define modified SMARTS patterns
    # Loosening the complete glucose requirement to a partial glucose or typical hexose structure
    partial_glucose_pattern = Chem.MolFromSmarts("[C@H]1([C@H]([C@H](O)[C@H](O)[C@H](O1)O)CO)S")
    # Pattern for sulfur-linkage with possible variation
    sulfur_linkage_pattern = Chem.MolFromSmarts("SC[C@H]1OC[C@@H](O)[C@@H](O)C1")
    # Adjusted pattern for N-sulfooxyimino group to ensure flexibility in its representation
    sulfooxyimino_group_pattern = Chem.MolFromSmarts("[N]=[O]S(=O)(=O)[O-]")

    # Check for essential glucosinolate components
    if not mol.HasSubstructMatch(partial_glucose_pattern):
        return False, "No recognizable glucose scaffold found"
    
    if not mol.HasSubstructMatch(sulfur_linkage_pattern):
        return False, "No identifiable sulfur linkage pattern found"

    if not mol.HasSubstructMatch(sulfooxyimino_group_pattern):
        return False, "Missing N-sulfooxyimino functional group"

    # If all key components are present, classify as glucosinolate
    return True, "Contains recognizable glucose scaffold, sulfur linkage, and N-sulfooxyimino functional group consistent with glucosinolates"