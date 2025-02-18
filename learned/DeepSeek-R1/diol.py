"""
Classifies: CHEBI:23824 diol
"""
"""
Classifies: CHEBI:23824 diol
"""
from rdkit import Chem

def is_diol(smiles: str):
    """
    Determines if a molecule is a diol based on its SMILES string.
    A diol contains exactly two hydroxy groups (-OH) regardless of their chemical environment.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find all hydroxyl groups (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    total_hydroxyl = len(mol.GetSubstructMatches(hydroxyl_pattern))
    
    # Subtract hydroxyls that are part of carboxylic acids
    carboxylic_acid_pattern = Chem.MolFromSmarts("[OH]-C(=O)")
    carb_oh = len(mol.GetSubstructMatches(carboxylic_acid_pattern))
    
    # Subtract hydroxyls that are part of enols (C=C-OH)
    enol_pattern = Chem.MolFromSmarts("C=C-O")
    enol_oh = len(mol.GetSubstructMatches(enol_pattern))
    
    # Adjusted hydroxyl count
    adjusted_oh = total_hydroxyl - carb_oh - enol_oh
    
    if adjusted_oh == 2:
        return True, f"Contains exactly two hydroxyl groups (adjusted count: {adjusted_oh})"
    else:
        return False, f"Adjusted hydroxyl count: {adjusted_oh}, need exactly 2"