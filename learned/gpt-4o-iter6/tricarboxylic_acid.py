"""
Classifies: CHEBI:27093 tricarboxylic acid
"""
from rdkit import Chem

def is_tricarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a tricarboxylic acid based on its SMILES string.
    A tricarboxylic acid is defined as having exactly three distinct carboxylic acid groups (-COOH).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a tricarboxylic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define carboxylic acid group using SMARTS
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    
    # Count the number of carboxylic acid groups, focusing on uniqueness of the carboxyl carbon
    distinct_carbons = set()
    for match in carboxylic_acid_matches:
        carbon_idx = match[0]
        distinct_carbons.add(carbon_idx)
    
    # Check the number of distinct carboxyl carbons
    if len(distinct_carbons) != 3:
        return False, f"Contains {len(distinct_carbons)} distinct carboxylic acid groups, expected exactly 3"

    # Verify that carboxylic acids are not part of large cycles or unexpected linkages
    for idx in distinct_carbons:
        atom = mol.GetAtomWithIdx(idx)
        # Check connectivity to ensure each -COOH group is distinct and separate
        neighbors = [n.GetIdx() for n in atom.GetNeighbors() if mol.GetAtomWithIdx(n.GetIdx()).GetSymbol() != 'C']
        if len(neighbors) != 2:  # Apart from the two different O atoms in -COOH
            return False, "Carboxylic acids not independently linked"

    return True, "Contains exactly three distinct and independent carboxylic acid groups"