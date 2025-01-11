"""
Classifies: CHEBI:23824 diol
"""
from rdkit import Chem

def is_diol(smiles: str):
    """
    Determines if a molecule is a diol based on its SMILES string.
    A diol should contain exactly two hydroxy groups that are part of an aliphatic structure.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a simple diol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for hydroxy group pattern (-OH) not explicitly tied to aromaticity but emphasizing aliphatic connections
    hydroxy_pattern = Chem.MolFromSmarts("[CX4;!$(C@C(O)C)][OX2H]")  # Carbon with attached hydroxyl not in aromatic systems
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    
    # Count the number of hydroxy groups
    num_hydroxy_groups = len(hydroxy_matches)

    # A diol should have exactly two -OH groups
    if num_hydroxy_groups != 2:
        return False, f"Contains {num_hydroxy_groups} hydroxy groups, need exactly 2"

    # Additional check to confirm hydroxyls are part of aliphatic structure by looking at adjacent atoms
    for match in hydroxy_matches:
        hydroxyl_idx = match[1]  # Index of the hydroxyl oxygen
        connected_carbons = [
            neighbor for neighbor in mol.GetAtomWithIdx(hydroxyl_idx).GetNeighbors()
            if neighbor.GetAtomicNum() == 6  # Ensure connected to carbon atoms
        ]

        # Ensure this carbon is part of aliphatic chain or non-aromatic cycle within context
        # Allow structures forming non-aromatic rings or simply linear
        for carbon in connected_carbons:
            if carbon.GetIsAromatic():
                return False, "Hydroxyl group is tied within aromatic structure, not suitable for simple diol"

    return True, "Contains exactly two hydroxy groups forming a suitable diol"