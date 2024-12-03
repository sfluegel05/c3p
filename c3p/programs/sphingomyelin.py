"""
Classifies: CHEBI:64583 sphingomyelin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_sphingomyelin(smiles: str):
    """
    Determines if a molecule is a sphingomyelin.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingomyelin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphorylcholine group
    phosphorylcholine = Chem.MolFromSmarts("COP([O-])(=O)OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphorylcholine):
        return False, "No phosphorylcholine group found"

    # Check for amide linkage with fatty acid
    amide = Chem.MolFromSmarts("NC(=O)")
    if not mol.HasSubstructMatch(amide):
        return False, "No amide linkage found"

    # Check for sphingoid base
    sphingoid_base = Chem.MolFromSmarts("[C@H](O)[C@H](O)[C@H](COP([O-])(=O)OCC[N+](C)(C)C)NC(=O)")
    if not mol.HasSubstructMatch(sphingoid_base):
        return False, "No sphingoid base found"

    return True, "Molecule is a sphingomyelin"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:64583',
                          'name': 'sphingomyelin',
                          'definition': 'Any of a class of phospholipids in '
                                        'which the amino group of a sphingoid '
                                        'base is in amide linkage with one of '
                                        'several fatty acids, while the '
                                        'terminal hydroxy group of the '
                                        'sphingoid base is esterified to '
                                        'phosphorylcholine.',
                          'parents': [   'CHEBI:140325',
                                         'CHEBI:35284',
                                         'CHEBI:35786',
                                         'CHEBI:36700']},
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
    'num_true_positives': 6,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 25,
    'precision': 1.0,
    'recall': 0.1935483870967742,
    'f1': 0.3243243243243243,
    'accuracy': None}