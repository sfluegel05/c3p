"""
Classifies: CHEBI:50263 2-hydroxydicarboxylic acid
"""
from rdkit import Chem

def is_2_hydroxydicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxydicarboxylic acid based on its SMILES string.
    A 2-hydroxydicarboxylic acid has two carboxylic acid groups and a hydroxy group
    on the carbon atom at position alpha to at least one of the carboxy groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 2-hydroxydicarboxylic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Pattern for carboxylic acid (-C(=O)O)
    carboxy_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxy_matches = mol.GetSubstructMatches(carboxy_pattern)
    
    # Ensure we have at least two distinct carboxylic acid groups
    if len(carboxy_matches) < 2:
        return False, f"Expected at least 2 carboxylic acid groups but found {len(carboxy_matches)}"
    
    # Generalized pattern for alpha-hydroxy: HO-C(-)-C(=O)O
    alpha_hydroxy_pattern = Chem.MolFromSmarts("O[C]C(=O)O")
    alpha_hydroxy_matches = mol.GetSubstructMatches(alpha_hydroxy_pattern)
    
    # Check that the hydroxy group is on alpha to a carboxylic acid group
    for match in alpha_hydroxy_matches:
        # Atom 0 is the hydroxy oxygen, Atom 1 is the alpha carbon, Atom 2 is connected to carboxylic group
        alpha_carbon = match[1]
        
        # Check if there exists another carboxylic group not necessarily directly connected but in the neighborhood
        for neighbor in mol.GetAtomWithIdx(alpha_carbon).GetNeighbors():
            if neighbor.HasSubstructMatch(carboxy_pattern) and neighbor.GetIdx() not in match:
                return True, "Contains two carboxylic acid groups and a hydroxy group on the alpha carbon"

    return False, "No hydroxy group found on the alpha carbon to a carboxylic acid group"