"""
Classifies: CHEBI:25627 octadecadienoic acid
"""
"""
Classifies: octadecadienoic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdchem

def is_octadecadienoic_acid(smiles: str):
    """
    Determines if a molecule is an octadecadienoic acid based on its SMILES string.
    An octadecadienoic acid is defined as 'Any straight-chain, C18 polyunsaturated fatty acid having two C=C double bonds.'

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octadecadienoic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group at the end
    carboxy_pattern = Chem.MolFromSmarts("C(=O)[O;H1,-]")  # Carboxylic acid functional group
    if not mol.HasSubstructMatch(carboxy_pattern):
        return False, "No carboxylic acid group found"

    # Get the longest carbon chain (should be unbranched)
    chains = mol.GetSubstructMatches(Chem.MolFromSmarts("[C;H2,H1](-*)-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*"))  # Matches linear carbon chains
    longest_chain = max(chains, key=lambda chain: len(chain)) if chains else None

    if longest_chain is None:
        return False, "No carbon chain found"
    elif len(longest_chain) != 18:
        return False, f"Main carbon chain has {len(longest_chain)} carbons, expected 18"

    # Check for branches (atoms connected to the chain carbons that are not chain members or hydrogens)
    for atom_idx in longest_chain:
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() not in longest_chain and neighbor.GetAtomicNum() != 1:
                return False, "Branching detected in the carbon chain"

    # Count the number of C=C double bonds in the chain
    double_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == rdchem.BondType.DOUBLE:
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            if begin_idx in longest_chain and end_idx in longest_chain:
                double_bonds += 1

    if double_bonds != 2:
        return False, f"Found {double_bonds} C=C double bonds in the main chain, expected 2"

    return True, "Molecule is a straight-chain C18 fatty acid with two C=C double bonds"

__metadata__ = {   'chemical_class': {   'id': None,
                              'name': 'octadecadienoic acid',
                              'definition': 'Any straight-chain, C18 polyunsaturated fatty acid having two C=C double bonds.',
                              'parents': None},
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
        'num_true_positives': None,
        'num_false_positives': None,
        'num_true_negatives': None,
        'num_false_negatives': None,
        'num_negatives': None,
        'precision': None,
        'recall': None,
        'f1': None,
        'accuracy': None}