"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
"""
Classifies: CHEBI:????? omega-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_omega_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an omega-hydroxy fatty acid based on its SMILES string.
    An omega-hydroxy fatty acid has a carboxyl group at one end and a hydroxyl group at the opposite end of the carbon chain.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an omega-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for carboxylic acid group (-COOH)
    carboxylic_acid = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group"
    # Get the carboxylic acid carbon
    cooh_matches = mol.GetSubstructMatches(carboxylic_acid)
    cooh_carbon = cooh_matches[0][0]  # First match, carboxylic acid carbon
    
    # Calculate distance matrix to find the farthest atom from cooh_carbon
    distance_matrix = Chem.GetDistanceMatrix(mol)
    distances = distance_matrix[cooh_carbon]
    terminal_atom = distances.argmax()
    
    # Check if terminal atom is a carbon
    terminal_atom_obj = mol.GetAtomWithIdx(terminal_atom)
    if terminal_atom_obj.GetAtomicNum() != 6:
        return False, "Terminal atom is not carbon"
    
    # Check for hydroxyl group on terminal carbon
    hydroxyl_found = False
    for neighbor in terminal_atom_obj.GetNeighbors():
        if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() >= 1:
            hydroxyl_found = True
            break
    if not hydroxyl_found:
        return False, "No hydroxyl group on terminal carbon"
    
    # Check that the main chain is a straight line (no branching)
    # Get the path from cooh_carbon to terminal_atom
    path = Chem.GetShortestPath(mol, cooh_carbon, terminal_atom)
    for atom_idx in path[1:-1]:  # Exclude first and last
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() != 6:
            return False, "Non-carbon atom in main chain"
        # Check for branching: each atom should have exactly two adjacent carbons in the path
        prev_idx = path[path.index(atom_idx) - 1]
        next_idx = path[path.index(atom_idx) + 1]
        neighbors_in_path = 0
        for n in atom.GetNeighbors():
            n_idx = n.GetIdx()
            if n_idx == prev_idx or n_idx == next_idx:
                neighbors_in_path += 1
        if neighbors_in_path != 2:
            return False, "Branching in main chain"
    
    return True, "Carboxylic acid at one end, hydroxyl at the opposite end of a straight carbon chain"