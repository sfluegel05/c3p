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

    # Define more specific patterns for common GAG substructure components
    # Amino sugar patterns: N-acetylglucosamine or N-acetylgalactosamine
    amino_sugar_pattern1 = Chem.MolFromSmarts("NC([C@H]1CO[C@H](O)[C@@H](O)[C@H]1O)C=O")
    amino_sugar_pattern2 = Chem.MolFromSmarts("NC([C@@H]1CO[C@H](O)[C@@H](O)[C@H]1O)C=O")
    # Uronic acid pattern: a sugar with a carboxylic acid
    uronic_acid_pattern = Chem.MolFromSmarts("[C@H]1(C(O)=O)O[C@H](O)[C@@H](O)[C@H](O)[C@H]1O")
    # Sulfate pattern
    sulfate_pattern = Chem.MolFromSmarts("OS(=O)(=O)O")

    # Check for the presence of amino sugars and uronic acid units
    amino_matches1 = mol.GetSubstructMatches(amino_sugar_pattern1)
    amino_matches2 = mol.GetSubstructMatches(amino_sugar_pattern2)
    uronic_matches = mol.GetSubstructMatches(uronic_acid_pattern)
    
    # Check for sulfate groups (optional)
    sulfate_matches = mol.GetSubstructMatches(sulfate_pattern)

    # Evaluate if the structure corresponds to a GAG
    if (len(amino_matches1) >= 1 or len(amino_matches2) >= 1) and len(uronic_matches) >= 1:
        reason = "Contains repeating units of amino sugars and uronic acids"
        if len(sulfate_matches) > 0:
            reason += " with sulfation"
        return True, reason
    
    return False, "Does not contain recognizable glycosaminoglycan features"