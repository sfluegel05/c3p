"""
Classifies: CHEBI:35915 sterol ester
"""
"""
Classifies: CHEBI:37664 sterol ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.
    A sterol ester is a steroid ester obtained by formal condensation of the carboxy group of any carboxylic acid with the 3-hydroxy group of a sterol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sterol backbone pattern (tetracyclic ring system with a 3-hydroxy group)
    sterol_pattern = Chem.MolFromSmarts("[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C")
    if not mol.HasSubstructMatch(sterol_pattern):
        return False, "No sterol backbone found"

    # Look for ester group (-O-C(=O)-) attached to the 3-hydroxy group
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) == 0:
        return False, "No ester group found"

    # Check if the ester group is attached to the 3-hydroxy group of the sterol
    ester_attached_to_sterol = False
    for match in ester_matches:
        ester_atom = mol.GetAtomWithIdx(match[0])
        for neighbor in ester_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon atom
                for neighbor_neighbor in neighbor.GetNeighbors():
                    if neighbor_neighbor.GetAtomicNum() == 8 and neighbor_neighbor.GetIdx() in [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetIdx() in [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetIdx() in [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8]]]:
                        ester_attached_to_sterol = True
                        break
                if ester_attached_to_sterol:
                    break
        if ester_attached_to_sterol:
            break

    if not ester_attached_to_sterol:
        return False, "Ester group not attached to the 3-hydroxy group of the sterol"

    return True, "Contains sterol backbone with ester linkage at the 3-hydroxy group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:37664',
                          'name': 'sterol ester',
                          'definition': 'A steroid ester obtained by formal condensation of the carboxy group of any carboxylic acid with the 3-hydroxy group of a sterol.',
                          'parents': ['CHEBI:37663', 'CHEBI:37665']},
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
    'num_true_positives': 150,
    'num_false_positives': 4,
    'num_true_negatives': 182407,
    'num_false_negatives': 23,
    'num_negatives': None,
    'precision': 0.974025974025974,
    'recall': 0.8670520231213873,
    'f1': 0.9174311926605504,
    'accuracy': 0.9998521228585199}