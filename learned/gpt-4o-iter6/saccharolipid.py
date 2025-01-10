"""
Classifies: CHEBI:166828 saccharolipid
"""
from rdkit import Chem

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.
    A saccharolipid should contain both a carbohydrate moiety and a lipid component.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a saccharolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string provided."
    
    # Broader carbohydrate pattern to identify various sugar rings
    carb_patterns = [
        Chem.MolFromSmarts("[OH]-C1(-[OH])O[C@H](O)[C@@H](O)[C@@H](O)[C@H]1"),        # Common sugar rings
        Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@@H](O)C[C@H](O)[C@@H]1O")]            # other sugar structures

    # Broad lipid pattern to identify generic long chains with ester or amide groups
    lipid_patterns = [
        Chem.MolFromSmarts("C(=O)[O,N][C@H1]-[C,C@](C)[C,C,H]"),                    # ester or amide groups with chains
        Chem.MolFromSmarts("C(=O)O")                                               # general ester group
    ]

    # Check for presence of a carbohydrate component
    carb_present = any(mol.HasSubstructMatch(carb_pattern) for carb_pattern in carb_patterns)
    if not carb_present:
        return False, "No adequate carbohydrate moiety found."

    # Check for presence of a lipid component
    lipid_present = any(mol.HasSubstructMatch(lipid_pattern) for lipid_pattern in lipid_patterns)
    if not lipid_present:
        return False, "No adequate long hydrocarbon chains (lipid component) found."

    return True, "Contains both carbohydrate and lipid components, indicating a saccharolipid."