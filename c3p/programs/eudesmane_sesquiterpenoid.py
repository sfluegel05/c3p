"""
Classifies: CHEBI:62508 eudesmane sesquiterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_eudesmane_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a eudesmane sesquiterpenoid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a eudesmane sesquiterpenoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for 15 carbon atoms (sesquiterpenoid)
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'C']
    if len(carbon_atoms) != 15:
        return False, "Molecule does not have 15 carbon atoms"

    # Check for eudesmane skeleton
    # Eudesmane skeleton: bicyclo[4.4.0]decane with methyl groups at positions 4 and 10
    try:
        eudesmane_skeleton = Chem.MolFromSmarts('CC1CCC2CCCCC2C1')
        if not mol.HasSubstructMatch(eudesmane_skeleton):
            return False, "Molecule does not have a eudesmane skeleton"
    except:
        return None, None

    return True, "Molecule is a eudesmane sesquiterpenoid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:62508',
                          'name': 'eudesmane sesquiterpenoid',
                          'definition': 'Any sesquiterpenoid having a '
                                        'eudesmane skeleton.',
                          'parents': ['CHEBI:26658']},
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
    'num_true_negatives': 13,
    'num_false_negatives': 7,
    'precision': 1.0,
    'recall': 0.46153846153846156,
    'f1': 0.631578947368421,
    'accuracy': None}