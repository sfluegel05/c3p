"""
Classifies: CHEBI:2440 acrovestone
"""
from rdkit import Chem

def is_acrovestone(smiles: str):
    """
    Determines if a molecule is an acrovestone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acrovestone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for polyphenol structure
    phenol_substructures = Chem.MolFromSmarts('c1ccccc1O')
    if not mol.HasSubstructMatch(phenol_substructures):
        return False, "No polyphenol structure found"

    # Check for glucoside linkage
    glucoside_linkage = Chem.MolFromSmarts('[C@H]1(O)[C@H](O)[C@@H](O)[C@H](O)[C@H](O1)')
    if not mol.HasSubstructMatch(glucoside_linkage):
        return False, "No glucoside linkage found"

    # Check for methoxy groups (OCH3)
    methoxy_group = Chem.MolFromSmarts('CO')
    if not mol.HasSubstructMatch(methoxy_group):
        return False, "No methoxy group found"

    return True, "Molecule is classified as acrovestone"

# Example usage
smiles_list = [
    "C1(=C(C2=C(C(C(=CO2)C3=CC=C(C(=C3)OC)OC)=O)C(=C1)O[C@H]4[C@H](C([C@@H](C(O4)CO)O)O)O[C@H]5[C@H](C([C@H](C(O5)C)O)O)O)C)OC",
    "O1[C@@H](OC2=CC=3OC[C@@H](CC3C=C2)C4=CC=C(O)C=C4)[C@H](O)[C@@H](O)[C@H](O)[C@H]1C(O)=O",
    "O1C(OC=2C=3OC=C(C(=O)C3C=CC2O)C4=CC=C(OC)C=C4)C(O)C(O)C(O)C1",
    "O1C(C(O)C(O)C(O)C1OC2=C(C(O)=C3C(OC=C(C3=O)C4=C(O)C=C(O)C=C4)=C2)CC=C(C)C)CO",
    "O1C(C(OC(=O)C)C(O)C(OC(=O)C)C1OC2=C(OC)C=C3C(OC=C(C3=O)C4=CC=C(O)C=C4)=C2)COC(=O)C",
    "O=C1C2=C(OC=C1C3=CC=C(O)C=C3)C=C(O[C@@H]4OC(C(OC)[C@H]([C@@H]4O)O)CO)C=C2O",
    "O1C(C(O)C(O)C(O)C1OC=2C=C3OC=C(C4=CC=C(OC5OCC(O)(C5O)CO)C=C4)C(=O)C3=C(O)C2)CO",
    "O=C1C2=C(OC=C1C3=CC=C(O[C@@H]4O[C@H]([C@H](O)[C@H]([C@H]4O)OC)C)C=C3)C=C(O[C@@H]5O[C@H]([C@H](O)[C@H]([C@H]5O)O)C)C=C2O",
    "O1C([C@@H](O)[C@H](O)C(O)[C@@H]1OC2=C(OC)C(O)=C3C(OC=C(C3=O)C4=CC(O)=C(OC)C=C4)=C2)CO",
    "O1C(C(O)C(O)C(OC(=O)/C=C/C2=CC=C(O)C=C2)C1OC=3C=C4OC=C(C(=O)C4=C(O)C3)C5=CC=C(O)C=C5)CO"
]

for smiles in smiles_list:
    result, reason = is_acrovestone(smiles)
    print(f"SMILES: {smiles} -> {result}, Reason: {reason}")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:2440',
                          'name': 'acrovestone',
                          'definition': 'A polyphenol that is isolated from '
                                        'Acronychia pedunculata and exhibits '
                                        'moderate antioxidant and '
                                        'antityrosinase activities.',
                          'parents': [   'CHEBI:22187',
                                         'CHEBI:26195',
                                         'CHEBI:35618',
                                         'CHEBI:78840']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 3,
    'success': True,
    'best': True,
    'error': '',
    'stdout': 'SMILES: '
              'C1(=C(C2=C(C(C(=CO2)C3=CC=C(C(=C3)OC)OC)=O)C(=C1)O[C@H]4[C@H](C([C@@H](C(O4)CO)O)O)O[C@H]5[C@H](C([C@H](C(O5)C)O)O)O)C)OC '
              '-> True, Reason: Molecule is classified as acrovestone\n'
              'SMILES: '
              'O1[C@@H](OC2=CC=3OC[C@@H](CC3C=C2)C4=CC=C(O)C=C4)[C@H](O)[C@@H](O)[C@H](O)[C@H]1C(O)=O '
              '-> True, Reason: Molecule is classified as acrovestone\n'
              'SMILES: '
              'O1C(OC=2C=3OC=C(C(=O)C3C=CC2O)C4=CC=C(OC)C=C4)C(O)C(O)C(O)C1 -> '
              'False, Reason: No glucoside linkage found\n'
              'SMILES: '
              'O1C(C(O)C(O)C(O)C1OC2=C(C(O)=C3C(OC=C(C3=O)C4=C(O)C=C(O)C=C4)=C2)CC=C(C)C)CO '
              '-> True, Reason: Molecule is classified as acrovestone\n'
              'SMILES: '
              'O1C(C(OC(=O)C)C(O)C(OC(=O)C)C1OC2=C(OC)C=C3C(OC=C(C3=O)C4=CC=C(O)C=C4)=C2)COC(=O)C '
              '-> True, Reason: Molecule is classified as acrovestone\n'
              'SMILES: '
              'O=C1C2=C(OC=C1C3=CC=C(O)C=C3)C=C(O[C@@H]4OC(C(OC)[C@H]([C@@H]4O)O)CO)C=C2O '
              '-> True, Reason: Molecule is classified as acrovestone\n'
              'SMILES: '
              'O1C(C(O)C(O)C(O)C1OC=2C=C3OC=C(C4=CC=C(OC5OCC(O)(C5O)CO)C=C4)C(=O)C3=C(O)C2)CO '
              '-> True, Reason: Molecule is classified as acrovestone\n'
              'SMILES: '
              'O=C1C2=C(OC=C1C3=CC=C(O[C@@H]4O[C@H]([C@H](O)[C@H]([C@H]4O)OC)C)C=C3)C=C(O[C@@H]5O[C@H]([C@H](O)[C@H]([C@H]5O)O)C)C=C2O '
              '-> True, Reason: Molecule is classified as acrovestone\n'
              'SMILES: '
              'O1C([C@@H](O)[C@H](O)C(O)[C@@H]1OC2=C(OC)C(O)=C3C(OC=C(C3=O)C4=CC(O)=C(OC)C=C4)=C2)CO '
              '-> True, Reason: Molecule is classified as acrovestone\n'
              'SMILES: '
              'O1C(C(O)C(O)C(OC(=O)/C=C/C2=CC=C(O)C=C2)C1OC=3C=C4OC=C(C(=O)C4=C(O)C3)C5=CC=C(O)C=C5)CO '
              '-> True, Reason: Molecule is classified as acrovestone\n',
    'num_true_positives': 9,
    'num_false_positives': 2,
    'num_true_negatives': 8,
    'num_false_negatives': 1,
    'precision': 0.8181818181818182,
    'recall': 0.9,
    'f1': 0.8571428571428572,
    'accuracy': None}