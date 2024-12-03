"""
Classifies: CHEBI:23849 diterpenoid
"""
from rdkit import Chem

def is_diterpenoid(smiles: str):
    """
    Determines if a molecule is a diterpenoid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diterpenoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Calculate the molecular formula
    formula = Chem.rdMolDescriptors.CalcMolFormula(mol)

    # Check if the molecular formula corresponds to a diterpenoid (C20 skeleton)
    if 'C20' not in formula:
        return False, "Molecular formula does not contain C20 skeleton"

    # Check if the molecule is a terpenoid (contains oxygen)
    if 'O' not in formula:
        return False, "Molecule does not contain oxygen, not a terpenoid"

    # Check for common diterpenoid features (e.g., rings, double bonds, functional groups)
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 1:
        return False, "No rings found, not a typical diterpenoid structure"

    # Check for specific functional groups that are common in diterpenoids
    functional_groups = ['C=O', 'C-O', 'C-OH']
    has_functional_group = any(mol.HasSubstructMatch(Chem.MolFromSmarts(fg)) for fg in functional_groups)
    if not has_functional_group:
        return False, "No common diterpenoid functional groups found"

    return True, "Molecule is classified as a diterpenoid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:23849',
                          'name': 'diterpenoid',
                          'definition': 'Any terpenoid derived from a '
                                        'diterpene. The term includes '
                                        'compounds in which the C20 skeleton '
                                        'of the parent diterpene has been '
                                        'rearranged or modified by the removal '
                                        'of one or more skeletal atoms '
                                        '(generally methyl groups).',
                          'parents': ['CHEBI:26873']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[15:19:39] SMILES Parse Error: syntax error while parsing: '
             'CO[C@]1(C)CC\\C=C(C)\\C[C@H]2OC(=O)C(CC[C@](O)(C(C)C)\\C=C\x01)=C2\n'
             '[15:19:39] SMILES Parse Error: Failed parsing SMILES '
             "'CO[C@]1(C)CC\\C=C(C)\\C[C@H]2OC(=O)C(CC[C@](O)(C(C)C)\\C=C\x01)=C2' "
             'for input: '
             "'CO[C@]1(C)CC\\C=C(C)\\C[C@H]2OC(=O)C(CC[C@](O)(C(C)C)\\C=C\x01)=C2'\n"
             '[15:19:39] SMILES Parse Error: syntax error while parsing: '
             '[H][C@@]12OC(=O)C(=C)[C@]1([H])CC\\C(C)=C\\CC\\C(=C/CC\\C(C)=C\x02)C(O)=O\n'
             '[15:19:39] SMILES Parse Error: Failed parsing SMILES '
             "'[H][C@@]12OC(=O)C(=C)[C@]1([H])CC\\C(C)=C\\CC\\C(=C/CC\\C(C)=C\x02)C(O)=O' "
             'for input: '
             "'[H][C@@]12OC(=O)C(=C)[C@]1([H])CC\\C(C)=C\\CC\\C(=C/CC\\C(C)=C\x02)C(O)=O'\n"
             '[15:19:39] SMILES Parse Error: syntax error while parsing: '
             'CO[C@]1(C)CC\\C=C(C)\\C[C@H]2OC(=O)C(CC[C@@H](\\C=C\x01)C(C)C)=C2\n'
             '[15:19:39] SMILES Parse Error: Failed parsing SMILES '
             "'CO[C@]1(C)CC\\C=C(C)\\C[C@H]2OC(=O)C(CC[C@@H](\\C=C\x01)C(C)C)=C2' "
             'for input: '
             "'CO[C@]1(C)CC\\C=C(C)\\C[C@H]2OC(=O)C(CC[C@@H](\\C=C\x01)C(C)C)=C2'\n"
             '[15:19:39] SMILES Parse Error: syntax error while parsing: '
             'CC(=C)C1=C/C=C(C)/CC\\C=C(C)\\C[C@H]2OC(=O)C(CC\x01)=C2\n'
             '[15:19:39] SMILES Parse Error: Failed parsing SMILES '
             "'CC(=C)C1=C/C=C(C)/CC\\C=C(C)\\C[C@H]2OC(=O)C(CC\x01)=C2' for "
             'input: '
             "'CC(=C)C1=C/C=C(C)/CC\\C=C(C)\\C[C@H]2OC(=O)C(CC\x01)=C2'\n",
    'stdout': '',
    'num_true_positives': 66,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 168,
    'precision': 1.0,
    'recall': 0.28205128205128205,
    'f1': 0.44000000000000006,
    'accuracy': None}