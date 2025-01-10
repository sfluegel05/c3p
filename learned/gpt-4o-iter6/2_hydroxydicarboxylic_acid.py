"""
Classifies: CHEBI:50263 2-hydroxydicarboxylic acid
"""
from rdkit import Chem

def is_2_hydroxydicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxydicarboxylic acid based on its SMILES string.
    A 2-hydroxydicarboxylic acid must have two carboxylic acid groups and a hydroxy group
    on the alpha carbon to at least one of the carboxylic acids.

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
    carboxy_pattern = Chem.MolFromSmarts("C(=O)[O;H1,-1]")
    carboxy_matches = mol.GetSubstructMatches(carboxy_pattern)
    
    # Ensure we have at least two distinct carboxylic acid groups
    if len(carboxy_matches) < 2:
        return False, f"Expected at least 2 carboxylic acid groups but found {len(carboxy_matches)}"

    # Refine pattern for alpha-hydroxy position: HO-C-C(=O)
    # Evaluating stereochemistry encoded examples more effectively
    alpha_hydroxy_patterns = [
        Chem.MolFromSmarts("O[C;H1,H2][C](C(=O)[O,H,-1])")
    ]

    for pattern in alpha_hydroxy_patterns:
        for match in mol.GetSubstructMatches(pattern):
            hydroxy_position = mol.GetAtomWithIdx(match[0])
            alpha_carbon_position = mol.GetAtomWithIdx(match[1])

            # Ensure the alpha carbon is directly connected to an additional carboxylic group
            for neighbor in alpha_carbon_position.GetNeighbors():
                if neighbor.GetSmarts() == '[C]' and neighbor.HasSubstructMatch(carboxy_pattern):
                    return True, "Contains two carboxylic acid groups and a hydroxy group on the alpha carbon"

    return False, "No hydroxy group found on the alpha carbon to a carboxylic acid group"