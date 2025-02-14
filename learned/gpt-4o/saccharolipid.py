"""
Classifies: CHEBI:166828 saccharolipid
"""
from rdkit import Chem

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.
    A saccharolipid contains a carbohydrate moiety linked to a lipid component.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a saccharolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define patterns for diverse carbohydrate moieties
    carbohydrate_patterns = [
        Chem.MolFromSmarts("C1OC(O)C(O)C(O)C(O)C1"),  # Hexopyranose
        Chem.MolFromSmarts("C1OC(O)C(O)C(O)C1"),      # Pentopyranose
        Chem.MolFromSmarts("[CX4](O)[CX4](O)[CX4](O)"),# Linear sugar-like
        Chem.MolFromSmarts("OC[C@H](O)[C@H](O)C=O"),   # Open sugars
    ]

    # Check for the presence of carbohydrate moieties
    has_carbohydrate = any(mol.HasSubstructMatch(pat) for pat in carbohydrate_patterns)
    if not has_carbohydrate:
        return False, "No recognized carbohydrate moiety pattern found"

    # Define patterns for lipid moieties
    lipid_patterns = [
        Chem.MolFromSmarts("C(=O)[O,N][CX4][CX4]"),    # Fatty acid ester/amide
        Chem.MolFromSmarts("C(=O)[O,N]C[C;R0]{6,}"),   # Long chain ester/amide
        Chem.MolFromSmarts("C(=O)O[C;R0]{5,}"),        # Ester linkage
    ]

    # Check for the presence of lipid moieties
    has_lipid = any(mol.HasSubstructMatch(pat) for pat in lipid_patterns)
    if not has_lipid:
        return False, "No recognized lipid moiety pattern found"

    # Verify direct linkage between carbohydrate and lipid moieties
    for carb_pat in carbohydrate_patterns:
        carb_matches = mol.GetSubstructMatches(carb_pat)
        for lipid_pat in lipid_patterns:
            lipid_matches = mol.GetSubstructMatches(lipid_pat)
            for carb_match in carb_matches:
                for lipid_match in lipid_matches:
                    if set(carb_match).intersection(lipid_match):
                        return True, "Carbohydrate and lipid components are chemically linked"

    return False, "Recognized moieties were not found to be directly linked"

# Note: Classification accuracy might still be improved with more specific patterns for the distinct variants of saccharolipids.