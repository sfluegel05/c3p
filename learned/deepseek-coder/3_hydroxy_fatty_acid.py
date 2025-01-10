"""
Classifies: CHEBI:59845 3-hydroxy fatty acid
"""
"""
Classifies: CHEBI:XXXXX 3-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acid based on its SMILES string.
    A 3-hydroxy fatty acid has a carboxylic acid group and a hydroxyl group at the 3-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for hydroxyl group (-OH) at the 3-position relative to the carboxylic acid
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4][CX4][CX4]([OH])")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group at the 3-position found"

    # Verify that the hydroxyl group is at the 3-position relative to the carboxylic acid
    # We need to ensure that the hydroxyl group is on the third carbon from the carboxylic acid
    # This can be done by checking the distance between the carboxylic acid carbon and the hydroxyl carbon
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)

    for ca_match in carboxylic_acid_matches:
        ca_carbon = ca_match[0]  # The carbon in the carboxylic acid group
        for oh_match in hydroxyl_matches:
            oh_carbon = oh_match[2]  # The carbon with the hydroxyl group
            # Calculate the shortest path between the carboxylic acid carbon and the hydroxyl carbon
            path = Chem.GetShortestPath(mol, ca_carbon, oh_carbon)
            if len(path) == 4:  # Path length of 4 means 3 bonds between the two carbons
                break
        else:
            continue
        break
    else:
        return False, "Hydroxyl group is not at the 3-position relative to the carboxylic acid"

    # Check for a long carbon chain (at least 6 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 6:
        return False, "Carbon chain too short for a fatty acid"

    return True, "Contains a carboxylic acid group and a hydroxyl group at the 3-position with a long carbon chain"