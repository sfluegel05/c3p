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
    carboxylic_pattern = Chem.MolFromSmarts("C(=O)(O)")
    carboxylic_matches = mol.GetSubstructMatches(carboxylic_pattern)
    
    if len(carboxylic_matches) < 2:
        return False, "Less than two carboxylic acid groups found"
    
    # Define the alpha-hydroxy pattern
    alpha_hydroxy_pattern = Chem.MolFromSmarts("O[CH](C(=O)O)")
    alpha_hydroxy_matches = mol.GetSubstructMatches(alpha_hydroxy_pattern)

    if len(alpha_hydroxy_matches) < 1:
        return False, "No hydroxy group on alpha carbon found"
    
    return True, "Contains two carboxylic acid groups and a hydroxy group on the alpha carbon"

# Examples to test the function
print(is_2_hydroxydicarboxylic_acid("CC(C(O)=O)C(C)(O)C(O)=O")) # 2,3-dimethylmalic acid