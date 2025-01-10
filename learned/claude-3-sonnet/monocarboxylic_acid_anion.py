"""
Classifies: CHEBI:35757 monocarboxylic acid anion
"""
"""
Classifies: monocarboxylic acid anion
A carboxylic acid anion formed when the carboxy group of a monocarboxylic acid is deprotonated.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a monocarboxylic acid anion based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a monocarboxylic acid anion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carboxylate group pattern (-C([O-])=O)
    carboxylate_pattern = Chem.MolFromSmarts('[C](=[OX1])[O-]')
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    
    # Count carboxylate groups
    num_carboxylate = len(carboxylate_matches)
    if num_carboxylate == 0:
        return False, "No carboxylate group found"
    
    # Look for protonated carboxylic acid groups (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts('[CX3](=O)[OX2H1]')
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    num_carboxylic_acid = len(carboxylic_acid_matches)
    
    # Total number of carboxylic groups (protonated + deprotonated)
    total_carboxy_groups = num_carboxylate + num_carboxylic_acid
    
    if total_carboxy_groups > 1:
        return False, f"Found {total_carboxy_groups} carboxylic/carboxylate groups, should have exactly one"
    
    if num_carboxylate != 1:
        return False, "Must have exactly one deprotonated carboxylate group"
        
    if num_carboxylic_acid > 0:
        return False, "Should not have any protonated carboxylic acid groups"
    
    # Additional check to ensure the carboxylate carbon is attached to something
    # (not just a free formate ion)
    for match in carboxylate_matches:
        carboxylate_carbon = mol.GetAtomWithIdx(match[0])
        neighbors = carboxylate_carbon.GetNeighbors()
        non_oxygen_neighbors = [n for n in neighbors if n.GetAtomicNum() != 8]
        if len(non_oxygen_neighbors) == 0:
            return False, "Carboxylate must be attached to a carbon chain/ring (not formate)"
    
    return True, "Contains exactly one deprotonated carboxylate group"