"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
from rdkit import Chem

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule is a glycosaminoglycan based on its SMILES string.
    Glycosaminoglycans must contain repeating disaccharide units and amino sugar features.

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
    
    # SMARTS patterns for sugar rings and amino substitutions
    sugar_ring_pattern = Chem.MolFromSmarts("C1([OH])[CH2][CH][OH][CH][OH][CH2]1")
    amino_sugar_pattern = Chem.MolFromSmarts("C1[CH][NH2][CH2][CH][OH][CH2]1")
    
    # Search for recognizable sugar and amino sugar patterns
    sugar_matches = mol.GetSubstructMatches(sugar_ring_pattern)
    amino_sugar_matches = mol.GetSubstructMatches(amino_sugar_pattern)

    # Require at least a certain number of amino sugars to classify as glycosaminoglycan
    if len(amino_sugar_matches) >= 3:
        return True, f"Contains {len(amino_sugar_matches)} amino sugar residues, suggesting glycosaminoglycan"

    return False, "Does not contain recognizable glycosaminoglycan features"