"""
Classifies: CHEBI:53666 kaurane diterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_kaurane_diterpenoid(smiles: str):
    """
    Determines if a molecule is a kaurane diterpenoid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a kaurane diterpenoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Kaurane skeleton is a specific tetracyclic structure. 
    # We will use SMARTS pattern to identify it.
    kaurane_smarts = 'C1CC2CCC3(C1)C(C4CCC(C4)C3C2)C'
    kaurane_pattern = Chem.MolFromSmarts(kaurane_smarts)
    
    if mol.HasSubstructMatch(kaurane_pattern):
        return True, "Molecule contains a kaurane skeleton"
    else:
        return False, "Molecule does not contain a kaurane skeleton"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:53666',
                          'name': 'kaurane diterpenoid',
                          'definition': 'A diterpenoid compound having a '
                                        'kaurane skeleton.',
                          'parents': ['CHEBI:23849']},
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
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 11,
    'num_false_negatives': 11,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}