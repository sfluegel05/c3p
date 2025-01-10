"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
"""
Classifies: CHEBI:27388 very long-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acid based on its SMILES string.
    A very long-chain fatty acid is defined as a fatty acid with a chain length greater than C22.
    Ultra-long-chain fatty acids have a chain length greater than C27.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a very long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add hydrogens to accurately find aliphatic chains
    mol = Chem.AddHs(mol)

    # Look for carboxylic acid group (-COOH)
    carboxylic_acid = Chem.MolFromSmarts('C(=O)[O;H1]')
    matches = mol.GetSubstructMatches(carboxylic_acid)
    if len(matches) == 0:
        return False, "No terminal carboxylic acid group found"

    # Identify all aliphatic chains ending with carboxylic acid
    max_chain_length = 0
    for match in matches:
        carboxyl_carbon_idx = match[0]
        carboxyl_oxygen_idx = match[2]

        # Check if carboxyl carbon is terminal (only one single bond)
        carboxyl_carbon = mol.GetAtomWithIdx(carboxyl_carbon_idx)
        single_bond_neighbors = [nbr for nbr in carboxyl_carbon.GetNeighbors() if mol.GetBondBetweenAtoms(carboxyl_carbon_idx, nbr.GetIdx()).GetBondType() == rdchem.BondType.SINGLE]
        if len(single_bond_neighbors) != 1:
            continue  # Not a terminal carboxyl group

        # Start from the carbon next to the carboxyl carbon
        start_atom_idx = single_bond_neighbors[0].GetIdx()

        # Find the longest linear carbon chain starting from start_atom_idx
        chain_atoms = []
        visited = set()

        def dfs(atom_idx):
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() != 6:
                return 0
            if atom_idx in visited:
                return 0
            visited.add(atom_idx)

            # Exclude atoms in rings
            if atom.IsInRing():
                return 0

            length = 1  # Current carbon
            max_length = 1
            for neighbor in atom.GetNeighbors():
                nbr_idx = neighbor.GetIdx()
                bond = mol.GetBondBetweenAtoms(atom_idx, nbr_idx)
                if bond.GetBondType() == rdchem.BondType.SINGLE and neighbor.GetAtomicNum() == 6:
                    length = 1 + dfs(nbr_idx)
                    if length > max_length:
                        max_length = length
            return max_length

        chain_length = dfs(start_atom_idx)
        if chain_length > max_chain_length:
            max_chain_length = chain_length

    if max_chain_length == 0:
        return False, "No suitable aliphatic chain found"

    total_chain_length = max_chain_length + 1  # Include carboxyl carbon

    # Check if chain length exceeds 22 carbons
    if total_chain_length > 22:
        if total_chain_length > 27:
            return True, f"Ultra-long-chain fatty acid with chain length C{total_chain_length}"
        else:
            return True, f"Very long-chain fatty acid with chain length C{total_chain_length}"
    else:
        return False, f"Chain length is C{total_chain_length}, not greater than C22"

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:27388',
                          'name': 'very long-chain fatty acid',
                          'definition': 'A fatty acid which has a chain length greater than C22. Very long-chain fatty acids which have a chain length greater than C27 are also known as ultra-long-chain fatty acids.',
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
    'accuracy': None}