"""
Classifies: CHEBI:166828 saccharolipid
"""
"""
Classifies: saccharolipid
"""
from rdkit import Chem
from rdkit.Chem import rdFMCS

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.
    A saccharolipid is a lipid where fatty acyl chains are directly attached to a sugar backbone via ester or amide linkages.

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
    
    # Define SMARTS patterns for various sugar moieties commonly found in saccharolipids
    sugar_smarts = [
        "[C@H]1(O)[C@@H](O)[C@H](O)[C@@H](CO)[C@H](O)[C@H]1O",  # glucose
        "[C@H]1(O)[C@@H](O)[C@H](O)[C@H](O[C@H]1CO)CO",          # glucosamine
        "O=C([C@@H](O)CO)[C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O",  # Kdo
    ]
    sugar_patterns = [Chem.MolFromSmarts(s) for s in sugar_smarts]
    
    # Search for sugar moieties
    has_sugar = False
    for pat in sugar_patterns:
        if mol.HasSubstructMatch(pat):
            has_sugar = True
            break
    if not has_sugar:
        return False, "No carbohydrate moiety (sugar ring) found"
    
    # Define SMARTS patterns for fatty acyl chains attached via ester linkage to sugar oxygen
    fatty_acyl_ester_to_sugar = Chem.MolFromSmarts("O[C;R][C](=O)[C;$(C([CH2])[CH2])]")  # Ester linkage to sugar oxygen
    # Define SMARTS patterns for fatty acyl chains attached via amide linkage to sugar nitrogen
    fatty_acyl_amide_to_sugar = Chem.MolFromSmarts("N[C;R][C](=O)[C;$(C([CH2])[CH2])]")  # Amide linkage to sugar nitrogen
    
    # Check for ester or amide linkages to sugar
    ester_matches = mol.GetSubstructMatches(fatty_acyl_ester_to_sugar)
    amide_matches = mol.GetSubstructMatches(fatty_acyl_amide_to_sugar)
    if not (ester_matches or amide_matches):
        return False, "No fatty acyl chains attached to sugar via ester or amide linkages found"

    # Exclude sphingolipids by checking for sphingosine backbone
    sphingosine_pattern = Chem.MolFromSmarts("C(CO)[NH][CH](O)CC=C")  # Simplified sphingosine backbone
    if mol.HasSubstructMatch(sphingosine_pattern):
        return False, "Contains sphingosine backbone, likely a glycosphingolipid"

    return True, "Contains sugar moiety with fatty acyl chains attached via ester or amide linkages"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'saccharolipid',
        'definition': 'Lipids that contain a carbohydrate moiety.',
    },
}