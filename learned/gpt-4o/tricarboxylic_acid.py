"""
Classifies: CHEBI:27093 tricarboxylic acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tricarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a tricarboxylic acid based on its SMILES string.
    A tricarboxylic acid contains exactly three distinct carboxy groups (-C(=O)OH).

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

    # Look for carboxylic acid pattern (C(=O)O[HYDROGEN])
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H]")
    carboxy_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)

    # Collect distinct carboxylic acids
    carboxy_atoms_sets = {frozenset(match) for match in carboxy_matches}

    # Count distinct carboxylic acid groups
    num_distinct_carboxy = len(carboxy_atoms_sets)

    if num_distinct_carboxy == 3:
        return True, "Contains exactly three distinct carboxylic acid groups"
    elif num_distinct_carboxy < 3:
        return False, f"Found {num_distinct_carboxy} distinct carboxylic acid groups, need exactly 3"
    else:
        return False, f"Found {num_distinct_carboxy} distinct carboxylic acid groups, need exactly 3"

# Example usage
smiles = "N[C@@H](CCC[C@H](NC(=O)CCC(O)=O)C(O)=O)C(O)=O"
result, reason = is_tricarboxylic_acid(smiles)
print(result, reason)