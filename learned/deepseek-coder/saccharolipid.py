"""
Classifies: CHEBI:166828 saccharolipid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.
    A saccharolipid is a lipid that contains a carbohydrate moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a saccharolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for lipid-like features (long hydrocarbon chains, ester bonds)
    # Look for ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 1:
        return False, "No ester groups found (lipid-like feature missing)"

    # Check for carbohydrate-like features (sugar rings, glycosidic bonds)
    # Look for sugar rings (pyranose or furanose)
    sugar_pattern = Chem.MolFromSmarts("[C@H]1[C@H](O)[C@H](O)[C@H](O)[C@H](O)O1")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if len(sugar_matches) < 1:
        return False, "No sugar rings found (carbohydrate-like feature missing)"

    # Check for glycosidic bonds (C-O-C between sugar and any other part)
    glycosidic_pattern = Chem.MolFromSmarts("[C@H]1[C@H](O)[C@H](O)[C@H](O)[C@H](O)O1.[OX2][CX4]")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
    if len(glycosidic_matches) < 1:
        return False, "No glycosidic bonds found (carbohydrate-like feature missing)"

    # Check for long aliphatic chains (lipid-like feature)
    aliphatic_pattern = Chem.MolFromSmarts("[CX4][CX4][CX4][CX4][CX4]")
    aliphatic_matches = mol.GetSubstructMatches(aliphatic_pattern)
    if len(aliphatic_matches) < 1:
        return False, "No long aliphatic chains found (lipid-like feature missing)"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 20:
        return False, "Too few carbons for saccharolipid"
    if o_count < 6:
        return False, "Too few oxygens for saccharolipid"

    return True, "Contains both lipid-like and carbohydrate-like features"