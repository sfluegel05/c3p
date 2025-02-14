"""
Classifies: CHEBI:2468 secondary alpha-hydroxy ketone
"""
"""
Classifies: secondary alpha-hydroxy ketone
"""
from rdkit import Chem

def is_secondary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is a secondary alpha-hydroxy ketone based on its SMILES string.
    A secondary alpha-hydroxy ketone is an alpha-hydroxy ketone in which the carbonyl group and 
    the hydroxy group are linked by a carbon bearing one hydrogen and one organyl group.

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

    # Define the SMARTS pattern for secondary alpha-hydroxy ketone
    pattern = Chem.MolFromSmarts('[O;H1]-[C;H1;D3](-[C;!H0])-[C;X3](=O)')
    matches = mol.GetSubstructMatches(pattern)
    
    if not matches:
        return False, "No secondary alpha-hydroxy ketone substructure found"

    # Further verify that the alpha carbon is secondary (connected to exactly two carbons and one hydrogen)
    for match in matches:
        oxygen_idx, alpha_c_idx, connected_c_idx, carbonyl_c_idx = match
        
        alpha_c = mol.GetAtomWithIdx(alpha_c_idx)
        if alpha_c.GetDegree() != 3:
            continue  # Not a secondary carbon (should be connected to three atoms)
        if alpha_c.GetTotalNumHs() != 1:
            continue  # Should have exactly one hydrogen

        # Check that alpha carbon is connected to hydroxyl oxygen, carbonyl carbon, and one carbon (organyl group)
        neighbor_atoms = alpha_c.GetNeighbors()
        elements = [atom.GetSymbol() for atom in neighbor_atoms]
        if elements.count('O') != 1:
            continue  # Should be connected to one hydroxyl oxygen
        if elements.count('C') != 2:
            continue  # Should be connected to two carbons (one carbonyl carbon and one organyl group)

        # Check that the connected carbon is not the carbonyl carbon
        connected_c = mol.GetAtomWithIdx(connected_c_idx)
        if connected_c.GetIdx() == carbonyl_c_idx:
            continue  # The connected carbon should not be the carbonyl carbon

        return True, "Secondary alpha-hydroxy ketone substructure found"

    return False, "No matching secondary alpha-hydroxy ketone substructure found"