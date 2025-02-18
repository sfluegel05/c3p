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
    A diol contains exactly two hydroxy groups (-OH) that are not part of carboxylic acids or phenols.

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
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    total_hydroxyl = len(hydroxyl_matches)
    
    # Find hydroxyls in carboxylic acids
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OH]")
    carboxylic_matches = len(mol.GetSubstructMatches(carboxylic_acid_pattern))
    
    # Find phenolic hydroxyl groups
    phenol_pattern = Chem.MolFromSmarts("[c]O[!H0]")
    phenol_matches = len(mol.GetSubstructMatches(phenol_pattern))
    
    # Calculate alcoholic hydroxyls
    alcoholic_hydroxyl = total_hydroxyl - carboxylic_matches - phenol_matches
    
    if alcoholic_hydroxyl == 2:
        return True, "Contains exactly two alcoholic hydroxyl groups"
    else:
        return False, f"Found {alcoholic_hydroxyl} alcoholic hydroxyl groups, need exactly 2"