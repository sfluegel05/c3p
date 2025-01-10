"""
Classifies: CHEBI:50699 oligosaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_oligosaccharide(smiles: str):
    """
    Determines if a molecule is an oligosaccharide based on its SMILES string.
    An oligosaccharide is a compound with monosaccharide units joined by glycosidic linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oligosaccharide, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns to search for hexose (6-membered sugar rings) units
    hexose_pattern = Chem.MolFromSmarts("[C@H]1(O)C(O)C(O)C(O)C(O)C1")
    
    # Search for at least two hexose units indicating oligosaccharide (can be adjusted for specificity)
    hexose_matches = mol.GetSubstructMatches(hexose_pattern)
    
    if len(hexose_matches) < 2:
        return False, f"Found {len(hexose_matches)} sugar units, need at least 2 for oligosaccharide"

    # Define glycosidic linkage pattern (O-link between sugars)
    glycosidic_linkage_pattern = Chem.MolFromSmarts("C-O-C")
    
    # Search for glycosidic linkages in the molecule
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_linkage_pattern)
    
    if len(glycosidic_matches) < 1:
        return False, "No glycosidic linkages found, important for oligosaccharide classification"
    
    return True, "Contains sufficient monosaccharide units linked by glycosidic bonds"