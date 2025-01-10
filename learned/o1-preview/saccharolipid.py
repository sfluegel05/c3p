"""
Classifies: CHEBI:166828 saccharolipid
"""
"""
Classifies: saccharolipid
"""
from rdkit import Chem

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
    
    # Define SMARTS patterns for sugar moieties (pyranose and furanose rings)
    sugar_patterns = [
        Chem.MolFromSmarts("C1OC(O)C(O)C(O)C1O"),  # Pyranose ring
        Chem.MolFromSmarts("C1OC(O)C(O)C1O"),      # Furanose ring
    ]
    
    # Search for sugar rings
    has_sugar = any(mol.HasSubstructMatch(pat) for pat in sugar_patterns)
    if not has_sugar:
        return False, "No carbohydrate moiety (sugar ring) found"
    
    # Define SMARTS pattern for fatty acyl chains attached via ester linkage to sugar
    # Fatty acyl chain: long aliphatic chain with carbonyl group
    # Ester linkage: [O][C](=O)[C]
    fatty_acyl_ester_pattern = Chem.MolFromSmarts("[C;!R][C;!R][C;!R][C;!R][C;!R][C;!R][C;!R][C;!R][C;!R][C;!R][C;!R][C;!R][CX3](=O)[O][#6]")
    
    # Define SMARTS pattern for fatty acyl chains attached via amide linkage to sugar
    # Amide linkage: [NX3][C](=O)[C]
    fatty_acyl_amide_pattern = Chem.MolFromSmarts("[C;!R][C;!R][C;!R][C;!R][C;!R][C;!R][C;!R][C;!R][C;!R][C;!R][C;!R][C;!R][CX3](=O)[NX3][#6]")
    
    # Check for fatty acyl chains connected to sugar via ester or amide linkages
    ester_bonds = mol.GetSubstructMatches(fatty_acyl_ester_pattern)
    amide_bonds = mol.GetSubstructMatches(fatty_acyl_amide_pattern)
    
    if not (ester_bonds or amide_bonds):
        return False, "No fatty acyl chains attached to sugar via ester or amide linkages found"
    
    return True, "Contains sugar moiety with fatty acyl chains attached via ester or amide linkages"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'saccharolipid',
        'definition': 'Lipids that contain a carbohydrate moiety.',
    },
}