"""
Classifies: CHEBI:3098 bile acid
"""
from rdkit import Chem

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid based on its SMILES string.
    Bile acids are hydroxy-5beta-cholanic acids with specific functional groups and stereochemistry.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for cholanic acid framework (basic steroid skeleton)
    steroid_skeleton = Chem.MolFromSmarts('[C@]12[C@]3([C@]4CC[C@@H]5[C@]1(CC[C@]5([C@@H]4C[C@H]3C2)C)C)C')
    if not mol.HasSubstructMatch(steroid_skeleton):
        return False, "Cholanic acid steroid skeleton not found"

    # Check for at least three hydroxyl groups
    hydroxy_pattern = Chem.MolFromSmarts('[OH]')
    if len(mol.GetSubstructMatches(hydroxy_pattern)) < 3:
        return False, "Insufficient number of hydroxy groups (less than three)"
    
    # Check for correct stereochemistry (5beta configuration)
    fivebeta_pattern = Chem.MolFromSmarts('[C@@H]1CCC[C@@H]2[C@]1(CCC3[C@@]2(CCC3)C)C')
    if not mol.HasSubstructMatch(fivebeta_pattern):
        return False, "5beta configuration not found"
    
    # Check for terminal carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts('C(=O)O')
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "Terminal carboxylic acid group not found"

    return True, "Valid bile acid: contains steroid core with hydroxy groups, 5beta stereochemistry, and terminal carboxylic group"

# Example test (Using one of the example SMILES)
# result, reason = is_bile_acid("O[C@@H]1[C@]2([C@]...")