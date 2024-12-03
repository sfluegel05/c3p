"""
Classifies: CHEBI:72823 glycerophosphoethanolamine zwitterion
"""
from rdkit import Chem
from rdkit.Chem import AllChem


def is_glycerophosphoethanolamine_zwitterion(smiles: str):
    """
    Determines if a molecule is a glycerophosphoethanolamine zwitterion.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophosphoethanolamine zwitterion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of the phosphate group
    phosphate_pattern = Chem.MolFromSmarts('P(=O)([O-])O')
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Check for the presence of the ethanolamine group
    ethanolamine_pattern = Chem.MolFromSmarts('NCCO')
    if not mol.HasSubstructMatch(ethanolamine_pattern):
        return False, "No ethanolamine group found"

    # Check for the presence of the glycerol backbone
    glycerol_pattern = Chem.MolFromSmarts('OCC(O)CO')
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Check for zwitterion characteristics: presence of both [NH3+] and [O-]
    zwitterion_pattern = Chem.MolFromSmarts('[NH3+].[O-]')
    if not mol.HasSubstructMatch(zwitterion_pattern):
        return False, "Molecule is not a zwitterion"

    return True, "Molecule is a glycerophosphoethanolamine zwitterion"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:72823',
                          'name': 'glycerophosphoethanolamine zwitterion',
                          'definition': 'A zwitterion obtained by transfer of '
                                        'a proton from the phosphate to the '
                                        'amino group of any '
                                        'glycerophosphoethanolamine.',
                          'parents': ['CHEBI:27369']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 29,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 3,
    'precision': 1.0,
    'recall': 0.90625,
    'f1': 0.9508196721311475,
    'accuracy': None}