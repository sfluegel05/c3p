"""
Classifies: CHEBI:36498 galactosylceramide
"""
from rdkit import Chem

def is_galactosylceramide(smiles: str):
    """
    Determines if a molecule is a galactosylceramide based on its SMILES string.
    Galactosylceramides are cerebrosides with a galactose head group and a sphingolipid backbone.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a galactosylceramide, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Correct the galactose portion detection
    galactose_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]1O")
    if galactose_pattern is None:
        return None, "Invalid SMARTS for galactose detection"
    if not mol.HasSubstructMatch(galactose_pattern):
        return False, "No recognizable galactose head group found"

    # Amide link pattern (C(=O)N)
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide linkage found"
    
    # Sphingolipid backbone detection
    # General pattern for sphingolipids includes [C@H](O)[C@H](COX), where X indicates the linkage to the galactose
    sphingolipid_pattern = Chem.MolFromSmarts("C(=O)N[C@@H](CO)C=C")
    if sphingolipid_pattern is None:
        return None, "Invalid SMARTS for sphingolipid detection"
    if not mol.HasSubstructMatch(sphingolipid_pattern):
        return False, "No sphingolipid backbone found"
    
    # Checking elemental composition for validity of carbohydrates and lipid chain
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    if o_count < 6:
        return False, "Insufficient oxygen atoms for necessary sugar and amide functionalities"
    if c_count < 30:
        return False, "Insufficient carbon atoms for long lipid and sugar chains"

    return True, "Contains galactose head group, sphingolipid backbone, and necessary structural features"