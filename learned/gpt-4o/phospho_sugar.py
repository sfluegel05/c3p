"""
Classifies: CHEBI:33447 phospho sugar
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.
    A phospho sugar is a monosaccharide containing an alcoholic hydroxy group esterified with phosphoric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phospho sugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generalized sugar pattern: flexible to include common saccharides. This uses a more abstract representation:
    sugar_pattern = Chem.MolFromSmarts("[C&R](O)[C&R](O)[C&R](O)[C&R](O)[CR](=O,[CX4H1])[O&R]")

    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar ring structure found"
    
    # General phosphate ester linkage: broadened pattern matching
    phosphate_pattern = Chem.MolFromSmarts("O[P](=O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate ester linkage found"
    
    return True, "Contains sugar ring structure with a phosphate ester linkage"