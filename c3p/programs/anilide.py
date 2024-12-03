"""
Classifies: CHEBI:13248 anilide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_anilide(smiles: str):
    """
    Determines if a molecule is an anilide (any aromatic amide obtained by acylation of aniline).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anilide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of aniline structure (benzene ring with NH2 group)
    aniline_pattern = Chem.MolFromSmarts('c1ccccc1N')
    if not mol.HasSubstructMatch(aniline_pattern):
        return False, "No aniline structure found"

    # Check for the presence of an amide bond (C=O attached to N)
    amide_pattern = Chem.MolFromSmarts('C(=O)N')
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide bond found"

    # Check if the amide nitrogen is directly attached to the benzene ring of aniline
    anilide_pattern = Chem.MolFromSmarts('c1ccccc1NC(=O)')
    if not mol.HasSubstructMatch(anilide_pattern):
        return False, "Amide nitrogen is not directly attached to the benzene ring of aniline"

    return True, "Molecule is an anilide"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:13248',
                          'name': 'anilide',
                          'definition': 'Any aromatic amide obtained by '
                                        'acylation of aniline.',
                          'parents': ['CHEBI:22712', 'CHEBI:62733']},
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
    'num_true_positives': 32,
    'num_false_positives': 5,
    'num_true_negatives': 15,
    'num_false_negatives': 0,
    'precision': 0.8648648648648649,
    'recall': 1.0,
    'f1': 0.927536231884058,
    'accuracy': None}