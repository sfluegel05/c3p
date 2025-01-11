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
    
    # Define sphingoid base pattern: long chain with hydroxyl and amino groups
    sphingoid_pattern = Chem.MolFromSmarts("O[C@H]([C@@H](O)CO)C")
    
    # Define amide linkage pattern
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    
    # Define general glycan (sugar) moiety pattern using a broad pattern
    # This is meant to identify various monosaccharides linked to oxygen
    glycan_pattern = Chem.MolFromSmarts("O[C@H]1O[C@@H](CO)[C@H](O)[C@@H](O)[C@H]1O")

    # Check for sphingoid base
    if not mol.HasSubstructMatch(sphingoid_pattern):
        return False, "No sphingoid base found"

    # Check for amide linkage
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide linkage found"

    # Check for glycan (sugar) moiety
    if not mol.HasSubstructMatch(glycan_pattern):
        return False, "No glycan (sugar) moiety found"

    return True, "Contains sphingoid base, amide linkage, and glycan moiety"