"""
Classifies: CHEBI:2468 secondary alpha-hydroxy ketone
"""
"""
Classifies: CHEBI:35654 secondary alpha-hydroxy ketone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_secondary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is a secondary alpha-hydroxy ketone based on its SMILES string.
    Secondary alpha-hydroxy ketones have a carbonyl group and a hydroxy group linked
    to the same carbon, which is also connected to an organyl group and a hydrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the alpha-hydroxy ketone pattern: C(=O)C(O)
    pattern = Chem.MolFromSmarts("[C&D3](=O)[C&D3](O)[!#1]")  # Exclude cases where C(=O)C(O) is in a ring
    matches = mol.GetSubstructMatches(pattern)
    
    if not matches:
        return False, "No secondary alpha-hydroxy ketone group found"
    
    for match in matches:
        atom_idx = match
        atom = mol.GetAtomWithIdx(atom_idx)
        
        # Check if the carbon is attached to an organyl group and a hydrogen
        organyl_count = sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetDegree() > 1)
        hydrogen_count = sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 1)
        
        if organyl_count < 1 or hydrogen_count != 1:
            continue  # Skip this match if conditions are not met
        
        # Check for exactly 1 carbonyl and 1 hydroxyl group
        carbonyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalDegree() == 1)
        hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalDegree() == 2)
        
        if carbonyl_count != 1 or hydroxyl_count != 1:
            continue  # Skip this match if conditions are not met
        
        return True, "Contains a secondary alpha-hydroxy ketone group"
    
    return False, "No valid secondary alpha-hydroxy ketone group found"