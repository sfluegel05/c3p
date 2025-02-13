"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_mucopolysaccharide(smiles: str):
    """
    Determines if a molecule is likely to be a mucopolysaccharide based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule has characteristics of a mucopolysaccharide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define cyclic uronic acid pattern (5/6-member ring with carboxylate)
    cyclic_uronic_acid_pattern = Chem.MolFromSmarts("[R]1[C;R][C;R][O;R][C;R][C;R](=O)O1")

    # Define cyclic glycosamine pattern (sugar with attached nitrogen)
    cyclic_glycosamine_pattern = Chem.MolFromSmarts("[R]1[C;R][C;R][C;R][N;R][C;R][O;R]1")

    # Define a more targeted sulfate ester pattern attached to carbohydrates
    sulfate_pattern = Chem.MolFromSmarts("[$([OX2]S(=O)(=O)[O][CX4;R])]")

    # Check for the presence of patterns
    uronic_acid_matches = mol.GetSubstructMatches(cyclic_uronic_acid_pattern)
    glycosamine_matches = mol.GetSubstructMatches(cyclic_glycosamine_pattern)
    sulfate_matches = mol.HasSubstructMatch(sulfate_pattern)

    if uronic_acid_matches and glycosamine_matches:
        if sulfate_matches:
            return True, "Characteristics of mucopolysaccharide, including sulfate esterification"
        else:
            return True, "Characteristics of mucopolysaccharide without detectable sulfate groups"
    
    return False, "Missing key features of a mucopolysaccharide"