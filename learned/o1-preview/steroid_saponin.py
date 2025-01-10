"""
Classifies: CHEBI:61655 steroid saponin
"""
"""
Classifies: CHEBI:26397 steroid saponin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_steroid_saponin(smiles: str):
    """
    Determines if a molecule is a steroid saponin based on its SMILES string.
    A steroid saponin is any saponin derived from a hydroxysteroid.
    This function checks for the presence of a steroid nucleus with hydroxyl groups,
    and sugar moieties attached via O-glycosidic bonds.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a steroid saponin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define steroid nucleus SMARTS pattern (tetracyclic fused ring system)
    steroid_pattern = Chem.MolFromSmarts('C1CCC2C1CCC3C2CCC4C3CCCC4')
    if steroid_pattern is None:
        return False, "Invalid steroid SMARTS pattern"

    # Check if steroid nucleus is present
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid nucleus found"

    # Get the atoms in the steroid nucleus
    steroid_matches = mol.GetSubstructMatch(steroid_pattern)
    steroid_atoms = set(steroid_matches)

    # Check for hydroxyl groups on the steroid nucleus carbons
    hydroxyl_on_steroid = False
    for atom in mol.GetAtoms():
        if atom.GetIdx() in steroid_atoms and atom.GetAtomicNum() == 6:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 1:
                    hydroxyl_on_steroid = True
                    break
            if hydroxyl_on_steroid:
                break

    if not hydroxyl_on_steroid:
        return False, "No hydroxyl groups on steroid nucleus found"

    # Check for sugar moieties attached via O-glycosidic bonds
    # Look for oxygen atoms connected to both steroid nucleus and a sugar ring
    glycosidic_bond = False
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            if atom1.GetAtomicNum() == 8 and atom2.GetAtomicNum() == 6:
                if (atom2.GetIdx() in steroid_atoms and atom1.GetDegree() == 2):
                    # Check if oxygen is connected to a sugar ring
                    for nb in atom1.GetNeighbors():
                        if nb.GetAtomicNum() == 6 and nb.GetIdx() not in steroid_atoms:
                            glycosidic_bond = True
                            break
            elif atom2.GetAtomicNum() == 8 and atom1.GetAtomicNum() == 6:
                if (atom1.GetIdx() in steroid_atoms and atom2.GetDegree() == 2):
                    # Check if oxygen is connected to a sugar ring
                    for nb in atom2.GetNeighbors():
                        if nb.GetAtomicNum() == 6 and nb.GetIdx() not in steroid_atoms:
                            glycosidic_bond = True
                            break
            if glycosidic_bond:
                break

    if not glycosidic_bond:
        return False, "No sugar moieties attached via O-glycosidic bonds found"

    return True, "Contains steroid nucleus with hydroxyl groups and sugar moieties attached via O-glycosidic bonds"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26397',
                              'name': 'steroid saponin',
                              'definition': 'Any saponin derived from a hydroxysteroid.',
                              'parents': ['CHEBI:26023', 'CHEBI:35188']},
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