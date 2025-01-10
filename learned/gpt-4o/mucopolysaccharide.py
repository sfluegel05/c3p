"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
from rdkit import Chem

def is_mucopolysaccharide(smiles: str):
    """
    Determines if a molecule is a mucopolysaccharide based on its SMILES string.
    Mucopolysaccharides are polysaccharides composed of repeating units of sugars which may include uronic acids 
    and glycosamines, often being partially esterified with sulfuric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mucopolysaccharide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generic sugar pattern: considering sugar ring with variations
    sugar_like_pattern = Chem.MolFromSmarts("C1OC([H,O])C([H,O])C(O)C1")
    if not mol.HasSubstructMatch(sugar_like_pattern):
        return False, "No sugar-like units found, unlikely to be a mucopolysaccharide"

    # Look for presence of any esterification pattern, including sulfate esters
    any_ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H0,R0]")  # General ester pattern
    sulfate_pattern = Chem.MolFromSmarts("O[S](=O)(=O)[O]")  # Specific sulfate ester
        
    if mol.HasSubstructMatch(sulfate_pattern):
        match_ester_msg = "sulfate ester groups"
    elif mol.HasSubstructMatch(any_ester_pattern):
        match_ester_msg = "generic ester groups"
    else:
        match_ester_msg = "no explicit ester groups"

    return True, f"Contains sugar-like units and {match_ester_msg}, indicating potential mucopolysaccharide"