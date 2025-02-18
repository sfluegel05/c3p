"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
"""
Classifies: CHEBI:????? omega-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import numpy as np

def is_omega_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an omega-hydroxy fatty acid based on its SMILES string.
    An omega-hydroxy fatty acid has a carboxyl group at one end and a hydroxyl group at the opposite end of a straight carbon chain.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an omega-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for exactly one carboxylic acid group (-COOH)
    carboxylic_acid = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    cooh_matches = mol.GetSubstructMatches(carboxylic_acid)
    if len(cooh_matches) != 1:
        return False, f"Expected exactly one carboxylic acid group, found {len(cooh_matches)}"
    cooh_carbon = cooh_matches[0][0]  # Carbon in the COOH group
    
    # Calculate distance matrix to find the farthest atom from COOH carbon
    distance_matrix = Chem.GetDistanceMatrix(mol)
    distances = distance_matrix[cooh_carbon]
    terminal_atom = int(np.argmax(distances))  # Convert numpy int to Python int
    
    # Terminal atom must be a carbon
    terminal_carbon = mol.GetAtomWithIdx(terminal_atom)
    if terminal_carbon.GetAtomicNum() != 6:
        return False, "Terminal atom is not carbon"
    
    # Check for hydroxyl group on terminal carbon
    hydroxyl_found = False
    for neighbor in terminal_carbon.GetNeighbors():
        if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() >= 1:
            hydroxyl_found = True
            break
    if not hydroxyl_found:
        return False, "No hydroxyl group on terminal carbon"
    
    # Check straight chain between COOH and terminal carbon
    path = Chem.GetShortestPath(mol, cooh_carbon, terminal_atom)
    if len(path) < 3:  # At least COOH-C-OH (minimum chain length 2 carbons?)
        return False, "Chain too short"
    
    for i, atom_idx in enumerate(path):
        atom = mol.GetAtomWithIdx(atom_idx)
        # All atoms in path must be carbons (except possibly O in COOH)
        if i == 0:
            continue  # COOH carbon is already checked
        if atom.GetAtomicNum() != 6:
            return False, f"Non-carbon atom in main chain at position {i}"
        
        # Check for branching in the main chain
        if i == 0 or i == len(path)-1:
            continue  # Endpoints can have more connections
        prev_idx = path[i-1]
        next_idx = path[i+1]
        neighbors_in_path = 0
        for n in atom.GetNeighbors():
            n_idx = n.GetIdx()
            if n_idx == prev_idx or n_idx == next_idx:
                neighbors_in_path += 1
        # Main chain carbons should have exactly 2 neighbors in the path (prev and next)
        if neighbors_in_path != 2:
            return False, f"Branching in main chain at carbon {atom_idx}"
    
    return True, "Carboxylic acid at one end, hydroxyl at the opposite end of a straight carbon chain"