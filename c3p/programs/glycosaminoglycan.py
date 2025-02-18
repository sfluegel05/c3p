"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
from rdkit import Chem

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule is a glycosaminoglycan based on its SMILES string.
    Glycosaminoglycans are characterized by containing repeat units of sugars 
    with substantial amino (NH2) modifications.
    
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
    
    # Consider generalized sugar-like patterns: simply a 5 or 6 member ring with multiple hydroxyls
    sugar_pattern = Chem.MolFromSmarts("C1([OH])[CH2][CH2][OH][CH][OH][CH1]1")
    
    # Consider an amino group in place of one of the hydroxyls, without rigid specificity
    aminosugar_pattern = Chem.MolFromSmarts("C1([OH])[CH2][CH2]([NH2])[CH][OH][CH1]1")
    
    # Check for these key substructures
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    aminosugar_matches = mol.GetSubstructMatches(aminosugar_pattern)

    # Assume it's a glycosaminoglycan if we detect amino-sugars in reasonable quantity
    if len(aminosugar_matches) >= 2:
        return True, "Contains a substantial proportion of aminomonosaccharide residues"

    return False, "Does not contain recognizable glycosaminoglycan features"