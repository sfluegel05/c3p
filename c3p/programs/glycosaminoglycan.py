"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
from rdkit import Chem

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule is a glycosaminoglycan based on its SMILES string.
    Glycosaminoglycans are characterized by repeating disaccharide units consisting of an amino sugar 
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
    aminomonosaccharide_pattern = Chem.MolFromSmarts("[C;X4][N]")  # Simple representation for amino sugars
    uronic_acid_pattern = Chem.MolFromSmarts("[C](=O)[O-]")  # Carboxylate functionality for uronic acids
    sulfate_pattern = Chem.MolFromSmarts("O=S(=O)([O-])[O-]")  # Sulfate group

    # Validate if patterns are correctly defined
    if aminomonosaccharide_pattern is None or uronic_acid_pattern is None:
        return None, "Error in defining substructure pattern"

    # Check for the presence of amino sugars and uronic acid units
    amino_matches = mol.GetSubstructMatches(aminomonosaccharide_pattern)
    uronic_matches = mol.GetSubstructMatches(uronic_acid_pattern)
    
    # Optional check for sulfate groups
    sulfate_matches = mol.GetSubstructMatches(sulfate_pattern) if sulfate_pattern is not None else []

    # Evaluate if the structure corresponds to a GAG
    if len(amino_matches) >= 1 and len(uronic_matches) >= 1:
        reason = "Contains repeating units of amino sugars and uronic acids"
        if len(sulfate_matches) > 0:
            reason += " with sulfation"
        return True, reason
    
    return False, "Does not contain recognizable glycosaminoglycan features"