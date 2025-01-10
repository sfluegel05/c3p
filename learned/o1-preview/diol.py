"""
Classifies: CHEBI:23824 diol
"""
"""
Classifies: CHEBI:23824 diol
"""
from rdkit import Chem

def is_diol(smiles: str):
    """
    Determines if a molecule is a diol based on its SMILES string.
    A diol is a compound that contains two hydroxy groups, generally assumed to be, but not necessarily, alcoholic.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Initialize count of valid hydroxyl groups
    hydroxyl_count = 0

    # Iterate over atoms to find hydroxyl groups
    for atom in mol.GetAtoms():
        # Check if atom is oxygen
        if atom.GetAtomicNum() == 8:
            # Check if oxygen is connected to hydrogen
            has_hydrogen = False
            has_carbon = False
            carbonyl_bond = False
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 1:
                    has_hydrogen = True
                elif neighbor.GetAtomicNum() == 6:
                    has_carbon = True
                    carbon = neighbor
            # If oxygen is bonded to hydrogen and carbon
            if has_hydrogen and has_carbon:
                # Check if carbon is part of a carbonyl group (C=O)
                for bond in carbon.GetBonds():
                    bonded_atom = bond.GetOtherAtom(carbon)
                    if bonded_atom.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        carbonyl_bond = True
                        break
                if not carbonyl_bond:
                    hydroxyl_count += 1

    if hydroxyl_count >= 2:
        return True, f"Molecule contains at least two hydroxyl groups attached to non-carbonyl carbons"
    else:
        return False, f"Molecule contains {hydroxyl_count} valid hydroxyl groups, diols must have at least two"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:23824',
                              'name': 'diol',
                              'definition': 'A compound that contains two hydroxy '
                                            'groups, generally assumed to be, but not '
                                            'necessarily, alcoholic. Aliphatic diols are '
                                            'also called glycols.',
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
        'attempt': 1,
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