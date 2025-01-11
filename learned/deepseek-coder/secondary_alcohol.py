"""
Classifies: CHEBI:35681 secondary alcohol
"""
"""
Classifies: CHEBI:17855 secondary alcohol
"""
from rdkit import Chem

def is_secondary_alcohol(smiles: str):
    """
    Determines if a molecule is a secondary alcohol based on its SMILES string.
    A secondary alcohol is a compound in which a hydroxy group, -OH, is attached to a saturated carbon atom which has two other carbon atoms attached to it.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more flexible substructure pattern for a secondary alcohol
    # The pattern is: [C](-[C])(-[C])-[OH]
    secondary_alcohol_pattern = Chem.MolFromSmarts("[C;!$(C=O)][C;!$(C=O)]([C;!$(C=O)])[OH]")
    
    # Find all matches of the pattern in the molecule
    matches = mol.GetSubstructMatches(secondary_alcohol_pattern)
    
    # Check if any of the matches are valid secondary alcohols
    for match in matches:
        # Get the carbon atom with the hydroxyl group
        carbon_idx = match[1]
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)
        
        # Ensure the carbon is saturated (sp3 hybridized)
        if carbon_atom.GetHybridization() != Chem.HybridizationType.SP3:
            continue
        
        # Ensure the carbon is attached to exactly two other carbons
        carbon_neighbors = [n for n in carbon_atom.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(carbon_neighbors) == 2:
            return True, "Contains a hydroxyl group attached to a carbon with two other carbon atoms"
    
    return False, "Does not contain a hydroxyl group attached to a carbon with two other carbon atoms"