"""
Classifies: CHEBI:76529 glycerophosphoglycerophosphoglycerol(2-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_glycerophosphoglycerophosphoglycerol_2__(smiles: str):
    """
    Determines if a molecule is a glycerophosphoglycerophosphoglycerol(2-).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophosphoglycerophosphoglycerol(2-), False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of glycerol phosphate moieties
    glycerol_phosphate_pattern = Chem.MolFromSmarts('OCC(O)COP(=O)([O-])O')
    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "No glycerol phosphate moiety found"

    # Check for two glycerol phosphate moieties linked to a central glycerol
    central_glycerol_pattern = Chem.MolFromSmarts('OCC(O)COP(=O)([O-])OCC(O)COP(=O)([O-])O')
    if not mol.HasSubstructMatch(central_glycerol_pattern):
        return False, "No central glycerol linking two glycerol phosphate moieties found"

    # Check for esterification to fatty acids
    ester_pattern = Chem.MolFromSmarts('C(=O)O')
    esters = mol.GetSubstructMatches(ester_pattern)
    if len(esters) < 2:
        return False, "Less than two esterified fatty acids found"

    return True, "Molecule is a glycerophosphoglycerophosphoglycerol(2-)"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:76529',
                          'name': 'glycerophosphoglycerophosphoglycerol(2-)',
                          'definition': 'An anionic phospholipid composed of '
                                        'two molecules of glycerol phosphate '
                                        'covalently linked to a molecule of '
                                        'glycerol, and in which each of the '
                                        'glycerol phosphate moieties may be '
                                        'esterified to one or two fatty acids.',
                          'parents': ['CHEBI:62643']},
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
    'num_true_positives': 10,
    'num_false_positives': 1,
    'num_true_negatives': 9,
    'num_false_negatives': 0,
    'precision': 0.9090909090909091,
    'recall': 1.0,
    'f1': 0.9523809523809523,
    'accuracy': None}