"""
Classifies: CHEBI:57643 1,2-diacyl-sn-glycero-3-phosphocholine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_1_2_diacyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is a 1,2-diacyl-sn-glycero-3-phosphocholine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1,2-diacyl-sn-glycero-3-phosphocholine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of the phosphocholine group
    phosphocholine_pattern = Chem.MolFromSmarts("P(OCC[N+](C)(C)C)(=O)[O-]")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "Missing phosphocholine group"

    # Check for the presence of the glycerol backbone with 1,2-diacyl substitution
    glycerol_pattern = Chem.MolFromSmarts("O[C@H](COP([O-])(=O)OCC[N+](C)(C)C)COC(=O)C")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Missing 1,2-diacyl-sn-glycero backbone"

    # Check for two acyl chains attached to the glycerol backbone
    acyl_chain_pattern = Chem.MolFromSmarts("OC(=O)C")
    acyl_chains = mol.GetSubstructMatches(acyl_chain_pattern)
    if len(acyl_chains) < 2:
        return False, "Less than two acyl chains detected"

    return True, "Valid 1,2-diacyl-sn-glycero-3-phosphocholine"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:57643',
                          'name': '1,2-diacyl-sn-glycero-3-phosphocholine',
                          'definition': 'The conjugate base of a '
                                        '1,2-diacyl-sn-glycero-3-phosphocholine '
                                        'compound formed by deprotonation of '
                                        'the phosphate OH group.',
                          'parents': [   'CHEBI:35284',
                                         'CHEBI:36313',
                                         'CHEBI:64482']},
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
    'num_true_positives': 76,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 28,
    'precision': 1.0,
    'recall': 0.7307692307692307,
    'f1': 0.8444444444444443,
    'accuracy': None}