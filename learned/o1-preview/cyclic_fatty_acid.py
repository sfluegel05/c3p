"""
Classifies: CHEBI:59238 cyclic fatty acid
"""
"""
Classifies: cyclic fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_cyclic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a cyclic fatty acid based on its SMILES string.
    A cyclic fatty acid is defined as any fatty acid containing anywhere in its structure a ring of atoms.

    Criteria:
    - Contains a carboxylic acid group (–COOH) or carboxylate anion (–COO⁻)
    - Has an aliphatic chain connected to the carboxyl group (chain length can vary)
    - Contains at least one ring (aliphatic or aromatic) anywhere in the molecule

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a cyclic fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group or carboxylate anion (-C(=O)OH or -C(=O)O-)
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)[OH]')
    carboxylate_pattern = Chem.MolFromSmarts('C(=O)[O-]')
    has_carboxylic_acid = mol.HasSubstructMatch(carboxylic_acid_pattern)
    has_carboxylate = mol.HasSubstructMatch(carboxylate_pattern)
    if not (has_carboxylic_acid or has_carboxylate):
        return False, "No carboxylic acid group found"

    # Check for ring structures in the molecule
    if not mol.GetRingInfo().NumRings():
        return False, "No ring structures found in the molecule"

    # Find the carboxyl carbon atoms
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    carboxyl_carbons = [match[0] for match in carboxylic_acid_matches + carboxylate_matches]

    # Check for aliphatic chain connected to the carboxyl carbon
    has_chain = False
    for carboxyl_carbon_idx in carboxyl_carbons:
        # Get the atom corresponding to the carboxyl carbon
        carboxyl_carbon_atom = mol.GetAtomWithIdx(carboxyl_carbon_idx)
        # Get neighbors excluding oxygens (focus on carbon chain)
        for neighbor in carboxyl_carbon_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon atom
                # Use BFS to find the carbon chain length
                visited = set()
                queue = [(neighbor, 1)]
                max_chain_length = 0
                while queue:
                    current_atom, length = queue.pop(0)
                    idx = current_atom.GetIdx()
                    if idx in visited:
                        continue
                    visited.add(idx)
                    if current_atom.GetAtomicNum() != 6:
                        continue
                    max_chain_length = max(max_chain_length, length)
                    for nbr in current_atom.GetNeighbors():
                        if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                            queue.append((nbr, length + 1))
                if max_chain_length >= 4:
                    has_chain = True
                    break
        if has_chain:
            break

    if not has_chain:
        return False, "No aliphatic chain connected to carboxylic acid group"

    return True, "Molecule meets criteria for cyclic fatty acid: contains carboxylic acid group, aliphatic chain, and ring structure"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'cyclic fatty acid',
        'definition': 'Any fatty acid containing anywhere in its structure a ring of atoms.',
        'parents': []
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': 50,
        'max_negative_in_prompt': 20,
        'max_instances_in_prompt': 100,
        'test_proportion': 0.1
    },
    'message': None,
    'attempt': 3,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': None,
    'num_false_positives': None,
    'num_true_negatives': None,
    'num_false_negatives': None,
    'num_negatives': None,
    'precision': None,
    'recall': None,
    'f1': None,
    'accuracy': None
}