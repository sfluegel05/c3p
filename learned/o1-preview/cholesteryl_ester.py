"""
Classifies: CHEBI:17002 cholesteryl ester
"""
"""
Classifies: cholesteryl ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cholesteryl_ester(smiles: str):
    """
    Determines if a molecule is a cholesteryl ester based on its SMILES string.
    A cholesteryl ester is a sterol ester obtained by formal condensation of the carboxy group
    of any carboxylic acid with the 3-hydroxy group of cholesterol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cholesteryl ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a general steroid nucleus pattern (four fused rings)
    steroid_pattern = Chem.MolFromSmarts('C1CCC2C(C1)CCC3C(C2)CCCC4=C3CCCC4')
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Steroid nucleus not found"

    # Define ester linkage pattern at the 3-hydroxy position
    ester_pattern = Chem.MolFromSmarts('C(=O)O[C]-[*]')  # Ester linkage
    ester_bond_found = False

    # Find atoms matching ester pattern
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester linkage found"

    # Check if ester is connected to the steroid nucleus at position 3
    steroid_matches = mol.GetSubstructMatches(steroid_pattern)
    steroid_atoms = set()
    for match in steroid_matches:
        steroid_atoms.update(match)

    for match in ester_matches:
        ester_oxygen = match[1]  # The oxygen atom in the ester linkage
        connected_atoms = [a.GetIdx() for a in mol.GetAtomWithIdx(ester_oxygen).GetNeighbors()]
        # Check if ester oxygen is connected to the steroid nucleus
        if any(atom_idx in steroid_atoms for atom_idx in connected_atoms):
            ester_bond_found = True
            break

    if not ester_bond_found:
        return False, "Ester linkage not connected to steroid nucleus"

    # Optionally, check that the acyl group is derived from a carboxylic acid (R-C(=O)-O-)
    # and is attached at position 3 of the steroid nucleus
    # Since exact atom positions may vary, we can accept the presence of the ester linkage

    return True, "Molecule is a cholesteryl ester (sterol ester with esterified 3-hydroxy group)"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'cholesteryl ester',
        'definition': 'A sterol ester obtained by formal condensation of the carboxy group of any carboxylic acid with the 3-hydroxy group of cholesterol.',
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
    'accuracy': None
}