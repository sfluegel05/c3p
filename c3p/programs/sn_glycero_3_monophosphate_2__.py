"""
Classifies: CHEBI:76530 sn-glycero-3-monophosphate(2-)
"""
from rdkit import Chem

def is_sn_glycero_3_monophosphate_2__(smiles: str):
    """
    Determines if a molecule is an sn-glycero-3-monophosphate(2-).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an sn-glycero-3-monophosphate(2-), False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core structure of sn-glycero-3-monophosphate(2-)
    core_smarts = "O[C@H](COP([O-])([O-])=O)CO"
    core = Chem.MolFromSmarts(core_smarts)
    if core is None:
        return None, None
    if not mol.HasSubstructMatch(core):
        return False, "Molecule does not contain sn-glycero-3-monophosphate(2-) core"

    # Check for acyl, alkyl, or alkenyl groups at sn-1 and sn-2 positions
    sn1_acyl_pattern = Chem.MolFromSmarts("C(=O)O[C@H](COP([O-])([O-])=O)CO")
    sn2_acyl_pattern = Chem.MolFromSmarts("O[C@H](COC(=O)C)COP([O-])([O-])=O")
    sn1_alkyl_pattern = Chem.MolFromSmarts("CO[C@H](COP([O-])([O-])=O)CO")
    sn2_alkyl_pattern = Chem.MolFromSmarts("O[C@H](COCC)COP([O-])([O-])=O")
    sn1_alkenyl_pattern = Chem.MolFromSmarts("C=C[C@H](COP([O-])([O-])=O)CO")
    sn2_alkenyl_pattern = Chem.MolFromSmarts("O[C@H](COC=C)COP([O-])([O-])=O")

    if sn1_acyl_pattern is None or sn2_acyl_pattern is None or sn1_alkyl_pattern is None or sn2_alkyl_pattern is None or sn1_alkenyl_pattern is None or sn2_alkenyl_pattern is None:
        return None, None

    if (mol.HasSubstructMatch(sn1_acyl_pattern) or mol.HasSubstructMatch(sn1_alkyl_pattern) or mol.HasSubstructMatch(sn1_alkenyl_pattern)) and \
       (mol.HasSubstructMatch(sn2_acyl_pattern) or mol.HasSubstructMatch(sn2_alkyl_pattern) or mol.HasSubstructMatch(sn2_alkenyl_pattern)):
        return True, "Molecule matches sn-glycero-3-monophosphate(2-) with acyl, alkyl, or alkenyl groups at sn-1 and sn-2 positions"
    
    return False, "Molecule does not match sn-glycero-3-monophosphate(2-) with required acyl, alkyl, or alkenyl groups at sn-1 and sn-2 positions"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:76530',
                          'name': 'sn-glycero-3-monophosphate(2-)',
                          'definition': 'An anionic phospholipid having a '
                                        'phosphate group at sn-3 position of '
                                        'the glycerol backbone, and with a '
                                        'combination of one or two acyl '
                                        'groups, alkyl groups, or alkenyl '
                                        'groups attached at the sn-1 and sn-2 '
                                        'positions through ester, ether or '
                                        'vinyl linkages respectively.',
                          'parents': ['CHEBI:62643', 'CHEBI:78185']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 4,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 5,
    'num_false_positives': 0,
    'num_true_negatives': 10,
    'num_false_negatives': 5,
    'precision': 1.0,
    'recall': 0.5,
    'f1': 0.6666666666666666,
    'accuracy': None}