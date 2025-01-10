"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
"""
Classifies: tetradecanoate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule is a tetradecanoate ester based on its SMILES string.
    A tetradecanoate ester is an ester formed from tetradecanoic acid (myristic acid) and an alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetradecanoate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define ester SMARTS pattern: carbonyl carbon (=O), single bond to oxygen, single bond to carbon
    ester_pattern = Chem.MolFromSmarts('[C;X3](=O)[O;X2][C]')
    if ester_pattern is None:
        return False, "Failed to create ester SMARTS pattern"

    # Find ester groups
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester groups found"

    # For each ester group, check if acyl chain is tetradecanoyl (14 carbons including carbonyl carbon)
    for match in ester_matches:
        carbonyl_c_idx = match[0]  # Index of carbonyl carbon
        ester_o_idx = match[2]     # Index of ester oxygen

        # Traverse acyl chain starting from carbonyl carbon, excluding ester oxygen
        visited = set()
        atoms_to_visit = [carbonyl_c_idx]
        c_count = 0

        while atoms_to_visit:
            atom_idx = atoms_to_visit.pop()
            if atom_idx in visited:
                continue
            visited.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            atomic_num = atom.GetAtomicNum()
            if atomic_num == 6:
                c_count += 1  # Count carbon atoms
            elif atomic_num != 1:
                # If we encounter a non-carbon, non-hydrogen atom, this is not a simple acyl chain
                break
            # Add neighbors to visit, excluding ester oxygen and already visited atoms
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx == ester_o_idx:
                    continue  # Do not traverse towards ester oxygen
                if neighbor_idx not in visited:
                    atoms_to_visit.append(neighbor_idx)

        if c_count == 14:
            return True, "Contains tetradecanoate ester group"

    return False, "Does not contain tetradecanoate ester group"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'tetradecanoate ester',
        'definition': 'A fatty acid ester obtained by condensation of the carboxy group of tetradecanoic acid (also known as myristic acid) with a hydroxy group of an alcohol or phenol.',
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
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None
}