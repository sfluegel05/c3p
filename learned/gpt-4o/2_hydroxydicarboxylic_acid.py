"""
Classifies: CHEBI:50263 2-hydroxydicarboxylic acid
"""
from rdkit import Chem

def is_2_hydroxydicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxydicarboxylic acid based on its SMILES string.
    
    A 2-hydroxydicarboxylic acid contains two carboxylic acid groups and a hydroxy 
    group on the carbon atom at position alpha to the carboxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-hydroxydicarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the carboxylic acid pattern
    carboxylic_pattern = Chem.MolFromSmarts("C(=O)[O]")
    carboxylic_matches = mol.GetSubstructMatches(carboxylic_pattern)
    
    if len(carboxylic_matches) < 2:
        return False, "Found less than two carboxylic acid groups"
    
    # Identify potential alpha-carbon positions
    carboxylic_ids = [match[0] for match in carboxylic_matches]  # Carbon atom index in carboxylic acids
    
    # Check for alpha-hydroxy pattern: Carbon with -OH directly attached and connected to a carboxylic carbon
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':  # Only interested in carbon atoms
            neighbors = atom.GetNeighbors()
            # Check if the carbon has both a hydroxy group and a carboxy group neighbors
            has_oh = any(nb.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(atom.GetIdx(), nb.GetIdx()).GetBondTypeAsDouble() == 1.0 for nb in neighbors)
            has_cox = any(nb.GetIdx() in carboxylic_ids for nb in neighbors)
            if has_oh and has_cox:
                return True, "Contains two carboxylic acid groups and a hydroxy group on the alpha carbon"
    
    return False, "No hydroxy group on alpha carbon found"

# Examples to test the function
print(is_2_hydroxydicarboxylic_acid("CC(C(O)=O)C(C)(O)C(O)=O")) # 2,3-dimethylmalic acid