"""
Classifies: CHEBI:16366 anthocyanidin cation
"""
from rdkit import Chem

def is_anthocyanidin_cation(smiles: str):
    """
    Determines if a molecule is an anthocyanidin cation.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anthocyanidin cation, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of the flavylium core structure
    flavylium_pattern = Chem.MolFromSmarts('c1cc2c(cc1)C=[O+]c3ccccc3O2')
    if not mol.HasSubstructMatch(flavylium_pattern):
        return False, "No flavylium core structure found"
    
    # Check for oxygenated derivatives
    oxygen_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'O']
    if len(oxygen_atoms) < 2:
        return False, "Not enough oxygen atoms for oxygenated derivatives"
    
    # Check if the molecule is a cation
    if not any(atom.GetFormalCharge() > 0 for atom in mol.GetAtoms()):
        return False, "Molecule is not a cation"
    
    return True, "Molecule is an anthocyanidin cation"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:16366',
                          'name': 'anthocyanidin cation',
                          'definition': 'Any organic cation that is an aglycon '
                                        'of anthocyanin cation; they are '
                                        'oxygenated derivatives of flavylium '
                                        '(2-phenylchromenylium).',
                          'parents': ['CHEBI:25697', 'CHEBI:47916']},
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
    'num_true_negatives': 18,
    'num_false_negatives': 18,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}