"""
Classifies: CHEBI:35746 fatty aldehyde
"""
"""
Classifies: fatty aldehyde
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a fatty aldehyde based on its SMILES string.
    A fatty aldehyde is an aldehyde arising from reduction of the carboxylic acid group of a fatty acid,
    having a carbonyl group at one end of the carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty aldehyde, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify aldehyde groups (-CHO)
    aldehyde_pattern = Chem.MolFromSmarts("[C;H1](=O)[#6]")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)

    if not aldehyde_matches:
        return False, "No aldehyde group found"
    
    # Check each aldehyde group
    for match in aldehyde_matches:
        aldehyde_carbon_idx = match[0]
        aldehyde_carbon = mol.GetAtomWithIdx(aldehyde_carbon_idx)
        
        # Get neighbors of aldehyde carbon (should be one single bond to carbon)
        neighbors = aldehyde_carbon.GetNeighbors()
        for neighbor in neighbors:
            if neighbor.GetAtomicNum() == 6:  # Carbon
                # Check if the chain connected to this carbon is aliphatic
                chain_atom = neighbor
                # Use BFS to traverse the chain and check if all atoms are aliphatic (exclude aromatic)
                visited = set()
                atoms_queue = [chain_atom.GetIdx()]
                is_aliphatic_chain = True
                while atoms_queue:
                    atom_idx = atoms_queue.pop(0)
                    if atom_idx in visited:
                        continue
                    visited.add(atom_idx)
                    atom = mol.GetAtomWithIdx(atom_idx)
                    if atom.GetIsAromatic():
                        is_aliphatic_chain = False
                        break
                    for nbr in atom.GetNeighbors():
                        nbr_idx = nbr.GetIdx()
                        if nbr_idx not in visited and nbr.GetAtomicNum() > 1:
                            atoms_queue.append(nbr_idx)
                if not is_aliphatic_chain:
                    continue  # Not an aliphatic chain, check next aldehyde
                else:
                    # Check if the aldehyde group is at the end of the chain (terminal)
                    if aldehyde_carbon.GetDegree() > 1:
                        continue  # Aldehyde carbon has more than one neighbor, not terminal
                    else:
                        return True, "Contains terminal aldehyde group attached to an aliphatic chain"
    return False, "No terminal aldehyde group attached to an aliphatic chain found"


__metadata__ = {   'chemical_class': {   'id': None,
                              'name': 'fatty aldehyde',
                              'definition': 'An aldehyde formally arising from reduction of the '
                                            'carboxylic acid group of its corresponding fatty acid, '
                                            'having a carbonyl group at one end of the carbon chain.',
                              'parents': []},
        'config': {   'llm_model_name': 'lbl/claude-sonnet',
                      'f1_threshold': 0.8,
                      'max_attempts': 5,
                      'max_positive_instances': None,
                      'max_positive_to_test': None,
                      'max_negative_to_test': None,
                      'max_positive_in_prompt': 50,
                      'max_negative_in_prompt': 20,
                      'max_instances_in_prompt': 100,
                      'test_proportion': 0.1},
        'message': None,
        'attempt': 0,
        'success': True,
        'best': True,
        'error': '',
        'stdout': None,
        'num_true_positives': 120,
        'num_false_positives': 10,
        'num_true_negatives': 182400,
        'num_false_negatives': 20,
        'num_negatives': None,
        'precision': 0.9230769230769231,
        'recall': 0.8571428571428571,
        'f1': 0.888888888888889,
        'accuracy': 0.9998345678901234}