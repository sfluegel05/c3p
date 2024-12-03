"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_nucleoside_5__phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside 5'-phosphate.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside 5'-phosphate, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a ribose or deoxyribose sugar
    ribose = Chem.MolFromSmarts('O[C@H]1[C@H](O)[C@@H](O[C@@H]1COP(O)(O)=O)')
    deoxyribose = Chem.MolFromSmarts('O[C@H]1[C@H](O)[C@@H](CO[C@@H]1COP(O)(O)=O)')
    if not mol.HasSubstructMatch(ribose) and not mol.HasSubstructMatch(deoxyribose):
        return False, "No ribose or deoxyribose sugar found"

    # Check for the presence of a pyrimidine or purine base
    pyrimidine = Chem.MolFromSmarts('c1cncnc1')
    purine = Chem.MolFromSmarts('c1ncnc2ncnc12')
    if not mol.HasSubstructMatch(pyrimidine) and not mol.HasSubstructMatch(purine):
        return False, "No pyrimidine or purine base found"

    # Check for the presence of a phosphate group at the 5' position
    phosphate = Chem.MolFromSmarts('COP(O)(O)=O')
    if not mol.HasSubstructMatch(phosphate):
        return False, "No phosphate group at the 5' position found"

    return True, "Molecule is a nucleoside 5'-phosphate"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:16701',
                          'name': "nucleoside 5'-phosphate",
                          'definition': 'A ribosyl or deoxyribosyl derivative '
                                        'of a pyrimidine or purine base in '
                                        'which C-5 of the ribose ring is '
                                        'mono-, di-, tri- or '
                                        'tetra-phosphorylated.',
                          'parents': ['CHEBI:29075']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[19:46:51] WARNING: not removing hydrogen atom without '
             'neighbors\n',
    'stdout': '',
    'num_true_positives': 15,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 5,
    'precision': 1.0,
    'recall': 0.75,
    'f1': 0.8571428571428571,
    'accuracy': None}