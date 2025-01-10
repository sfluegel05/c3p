"""
Classifies: CHEBI:24402 glycosphingolipid
"""
from rdkit import Chem

def is_glycosphingolipid(smiles: str):
    """
    Determines if a molecule is a glycosphingolipid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycosphingolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define patterns
    sphingoid_pattern = Chem.MolFromSmarts("C[C@H](O)CC")
    amide_pattern = Chem.MolFromSmarts("CN(C=O)C")
    sugar_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@H](O)[C@H](O)[C@H]([C@H]1O)O")
    
    # Check for sphingoid base
    if not mol.HasSubstructMatch(sphingoid_pattern):
        return False, "No sphingoid base found"
    
    # Check for amide linkage
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide linkage found"
    
    # Check for sugar moiety
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No glycan (sugar) moiety found"

    return True, "Contains sphingoid base, amide linkage, and glycan moiety"