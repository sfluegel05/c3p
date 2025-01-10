"""
Classifies: CHEBI:166828 saccharolipid
"""
from rdkit import Chem

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.
    A saccharolipid is characterized by the presence of both a lipid and a carbohydrate moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a saccharolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for comprehensive carbohydrate moiety pattern
    carbohydrate_patterns = [
        Chem.MolFromSmarts('[O;R]C(O)C(O)'),  # Represents carbohydrate rings with multiple hydroxyls
        Chem.MolFromSmarts('O[C@H]1O[C@H]')   # Cyclic ether structure typical of sugars
    ]
    carbohydrate_found = any(mol.HasSubstructMatch(p) for p in carbohydrate_patterns)
    
    if not carbohydrate_found:
        return False, "No carbohydrate moiety found"

    # Look for comprehensive lipid moiety pattern
    lipid_patterns = [
        Chem.MolFromSmarts('C(=O)OCCCCC'),   # Carboxylic esters with long chains
        Chem.MolFromSmarts('CCCCCCCCC(=O)')  # Long aliphatic chains ending in carbonyl
    ]
    lipid_found = any(mol.HasSubstructMatch(p) for p in lipid_patterns)
    
    if not lipid_found:
        return False, "No lipid moiety found"
    
    # Check for established linkage between carbohydrate and lipid
    linkage_patterns = [
        Chem.MolFromSmarts('[C@H]([O][C@H](C)C)([O][C][C](=O))')  # Glycosidic linkages to lipid esters
    ]
    linkage_found = any(mol.HasSubstructMatch(p) for p in linkage_patterns)
    
    if not linkage_found:
        return False, "No linkage between carbohydrate and lipid found"

    return True, "Contains linked carbohydrate and lipid moieties, consistent with a saccharolipid"