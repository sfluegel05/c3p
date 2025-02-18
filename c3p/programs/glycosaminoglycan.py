"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
from rdkit import Chem

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule is a glycosaminoglycan based on its SMILES string.
    A GAG is characterized by having repeating aminomonosaccharide units.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycosaminoglycan, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Defining a pattern for an amino sugar, common in GAGs
    aminomonosaccharide_pattern = Chem.MolFromSmarts("[OX2][CX4][NX3]")
    
    # Check if pattern exists in the molecule
    if mol.HasSubstructMatch(aminomonosaccharide_pattern):
        return True, "Contains aminomonosaccharide units characteristic of GAGs"
    
    return False, "Does not contain recognizable glycosaminoglycan features"