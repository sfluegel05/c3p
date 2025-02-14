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

    # Define a more generalized pattern for a carbohydrate moiety: 
    # Monosaccharides can vary widely, so include multiple recognizable features.
    carbohydrate_patterns = [
        Chem.MolFromSmarts("C1OC(O)C(O)C(O)C(O)C1"),  # Hexopyranose pattern
        Chem.MolFromSmarts("C1OC(O)C(O)C(O)C1"),      # Pentose pattern
        Chem.MolFromSmarts("[CX4](O)[CX4](O)[CX4](O)"),# Linear sugar-like structure
        Chem.MolFromSmarts("OC[C@H](O)[C@H](O)C=O"),   # Open form of glucose
    ]

    if not any(mol.HasSubstructMatch(pat) for pat in carbohydrate_patterns):
        return False, "No recognized carbohydrate moiety pattern found"

    # Define a pattern for a lipid moiety: longer hydrocarbon chain, ester or other typical lipid linkage
    lipid_patterns = [
        Chem.MolFromSmarts("C(=O)[O,N][CX4][CX4]"),    # Fatty acid ester or amide
        Chem.MolFromSmarts("C(=O)[O,N]C[C;R0]{6,}"),   # Long chain ester/amide
        Chem.MolFromSmarts("C(=O)O[C;R0]{5,}"),        # Ester linkage
    ]

    if not any(mol.HasSubstructMatch(pat) for pat in lipid_patterns):
        return False, "No recognized lipid moiety pattern found"

    # Check for plausible linkages between detected carbohydrate and lipid parts
    for carb_pat in carbohydrate_patterns:
        carb_matches = mol.GetSubstructMatches(carb_pat)
        for lipid_pat in lipid_patterns:
            lipid_matches = mol.GetSubstructMatches(lipid_pat)
            for carb_match in carb_matches:
                for lipid_match in lipid_matches:
                    if set(carb_match).intersection(lipid_match):
                        return True, "Carbohydrate and lipid components are chemically linked"

    return False, "Recognized moieties were not found to be directly linked"

# Note: Due to structural diversity of saccharolipids, further expansion might still be needed for uncommon structures.