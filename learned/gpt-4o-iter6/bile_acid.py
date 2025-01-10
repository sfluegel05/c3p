"""
Classifies: CHEBI:3098 bile acid
"""
from rdkit import Chem

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid based on its SMILES string.
    Bile acids are hydroxy-5beta-cholanic acids occurring as sodium salts of their amides with glycine or taurine.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Flexible 5beta-cholanic acid core pattern
    core_pattern = Chem.MolFromSmarts('[C@]12[C@@H]3CC[C@@H]4[C@@H](CCC[C@H]4[C@H]3CCC1)CC2')
    if not mol.HasSubstructMatch(core_pattern):
        return False, "5-beta-cholanic nucleus not identified"

    # Check for OH groups (minimum of 2 hydroxy groups)
    hydroxy_pattern = Chem.MolFromSmarts('[OH]')
    hydroxyl_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if len(hydroxyl_matches) < 2:
        return False, "Insufficient hydroxy groups (less than two)"

    # Check for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts('C(=O)O')
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "Carboxylic acid group not identified"

    return True, "Valid bile acid: contains 5beta-cholanic nucleus, multiple hydroxyl groups, and a terminal carboxylic group"

# Example usage:
# result, reason = is_bile_acid("O[C@@H]1[C@]2([C@]...")