"""
Classifies: CHEBI:33563 glycolipid
"""
from rdkit import Chem

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid based on its SMILES string.
    A glycolipid typically consists of a lipid moiety and a carbohydrate moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of long lipid-like carbon chains
    lipid_pattern = Chem.MolFromSmarts("C(=O)CCCCC")  # Simple pattern to find long chains (lipid tail indicator)
    if not mol.HasSubstructMatch(lipid_pattern):
        return False, "No lipid-like structure found"

    # Look for carbohydrate (sugar) moieties, extend the search with common monosaccharide substructures
    sugar_patterns = [
        "O[C@H]1[C@H](O)[C@@H](O)[C@H](CO)O[C@H]1O",  # Glucose-like pattern
        "O[C@H]1[C@H](O)[C@H](O)[C@@H](CO)O[C@@H]1",  # Galactose-like pattern
        "O[C@H](CO)[C@H]1O[C@@H](O)[C@@H](O)[C@H](O)[C@H]1O",  # Another sugar pattern variant
        "O[C@@H]1(CO)O[C@H](O[C@H]1O)C"  # Simplified furanose/pyranose ring patterns
    ]

    sugar_found = False
    for pattern in sugar_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            sugar_found = True
            break

    if not sugar_found:
        return False, "No carbohydrate moiety found"

    return True, "Contains both lipid and carbohydrate moieties consistent with glycolipids"