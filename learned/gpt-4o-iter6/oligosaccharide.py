"""
Classifies: CHEBI:50699 oligosaccharide
"""
from rdkit import Chem

def is_oligosaccharide(smiles: str):
    """
    Determines if a molecule is an oligosaccharide based on its SMILES string.
    An oligosaccharide is a compound with monosaccharide units joined by glycosidic linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an oligosaccharide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS patterns for different sugar backbones (hexoses and pentoses)
    sugar_patterns = [
        Chem.MolFromSmarts("[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O"),  # Example for hexoses
        Chem.MolFromSmarts("[C@H](O)[C@H](O)[C@@H](O)[C@H]1O"),  # Example for pentoses
        Chem.MolFromSmarts("[C@@H](O)[C@@H](O)[C@H](O)[C@H]1O")
        # Additional patterns for other common sugars can be added here
    ]
    
    # Search for sugar substructures
    sugar_matches = 0
    sugar_atoms = set()
    for pattern in sugar_patterns:
        matches = mol.GetSubstructMatches(pattern)
        sugar_matches += len(matches)
        for match in matches:
            sugar_atoms.update(match)
    
    if sugar_matches < 2:
        return False, f"Found {sugar_matches} sugar units, need at least 2 for an oligosaccharide"
    
    # Look for glycosidic linkages (O-links between two sugars)
    linkage_pattern = Chem.MolFromSmarts("C-O-C")
    
    # Ensure that there is at least one glycosidic linkage connecting distinct sugar units
    glycosidic_match = any(mol.HasSubstructMatch(linkage_pattern))
    
    if not glycosidic_match:
        return False, "No glycosidic linkages found, important for oligosaccharide classification"
    
    return True, "Contains sufficient monosaccharide units linked by glycosidic bonds"