"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
from rdkit import Chem

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule is a glycosaminoglycan based on its SMILES string.
    A GAG is characterized by having repeating disaccharide units consisting of an amino sugar 
    and a uronic acid, often with sulfate groups.
    
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

    # Define the substructure patterns typical in GAGs
    aminomonosaccharide_pattern = Chem.MolFromSmarts("[#6][OX2H][#6]([NX3+0])")  # Amino sugar with potential acetylation
    uronic_acid_pattern = Chem.MolFromSmarts("[#6](=O)[#8]-")  # Carboxylate functionality
    sulfate_pattern = Chem.MolFromSmarts("[$([OX2]S(=O)(=O)[O-])]")  # Sulfate group

    # Check for the presence of a repeating unit using substructure search
    amino_matches = mol.GetSubstructMatches(aminomonosaccharide_pattern)
    uronic_matches = mol.GetSubstructMatches(uronic_acid_pattern)
    sulfate_matches = mol.GetSubstructMatches(sulfate_pattern)

    # Evaluate if the structure corresponds to a GAG
    if len(amino_matches) > 1 and len(uronic_matches) > 1 and len(sulfate_matches) > 0:
        return True, "Contains repeating units of amino sugars and uronic acids with sulfation"
    
    return False, "Does not contain recognizable glycosaminoglycan features"