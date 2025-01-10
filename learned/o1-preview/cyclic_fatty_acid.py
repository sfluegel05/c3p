"""
Classifies: CHEBI:59238 cyclic fatty acid
"""
"""
Classifies: cyclic fatty acid
"""
from rdkit import Chem

def is_cyclic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a cyclic fatty acid based on its SMILES string.
    A cyclic fatty acid is defined as a long-chain aliphatic carboxylic acid containing anywhere in its structure a ring of atoms.

    Criteria:
    - Contains a carboxylic acid group (–COOH) or carboxylate anion (–COO⁻)
    - Has a long aliphatic chain (typically at least 8 carbons) connected to the carboxyl group
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

    # Function to count the length of the aliphatic carbon chain connected to the carboxyl carbon
    def count_aliphatic_chain(atom_idx, visited):
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom_idx in visited:
            return 0
        visited.add(atom_idx)
        if atom.GetAtomicNum() != 6 or atom.GetIsAromatic():
            return 0
        max_length = 1
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            bond = mol.GetBondBetweenAtoms(atom_idx, neighbor_idx)
            # Exclude non-single bonds and bonds that lead back to the carboxyl group
            if bond.GetBondType() != Chem.BondType.SINGLE:
                continue
            if neighbor.GetAtomicNum() != 6:
                continue
            if neighbor_idx in carboxyl_carbons:
                continue
            chain_length = 1 + count_aliphatic_chain(neighbor_idx, visited.copy())
            if chain_length > max_length:
                max_length = chain_length
        return max_length

    # Check for a long aliphatic chain connected to the carboxyl carbon
    has_long_chain = False
    for carboxyl_carbon_idx in carboxyl_carbons:
        # Skip atoms that are part of a ring
        if mol.GetAtomWithIdx(carboxyl_carbon_idx).IsInRing():
            continue
        chain_length = 0
        for neighbor in mol.GetAtomWithIdx(carboxyl_carbon_idx).GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor.GetAtomicNum() == 6 and not neighbor.GetIsAromatic():
                length = count_aliphatic_chain(neighbor_idx, set([carboxyl_carbon_idx]))
                if length > chain_length:
                    chain_length = length
        if chain_length >= 8:
            has_long_chain = True
            break

    if not has_long_chain:
        return False, "No sufficiently long aliphatic chain connected to carboxylic acid group"

    return True, "Molecule meets criteria for cyclic fatty acid: carboxylic acid group with long aliphatic chain and ring structure present"

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
    'attempt': 2,
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