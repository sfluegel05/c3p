"""
Classifies: CHEBI:83264 beta-D-glucosylceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_beta_D_glucosylceramide(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosylceramide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-glucosylceramide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of the beta-D-glucosyl part
    beta_D_glucosyl_smiles = "O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O1"
    beta_D_glucosyl = Chem.MolFromSmiles(beta_D_glucosyl_smiles)
    if beta_D_glucosyl is None:
        return False, "Invalid beta-D-glucosyl SMILES string"
    
    if not mol.HasSubstructMatch(beta_D_glucosyl):
        return False, "No beta-D-glucosyl group found"

    # Check for the presence of the ceramide part
    ceramide_smiles = "NC(CO)C(O)CCCCCCCCCC(C)C"
    ceramide = Chem.MolFromSmiles(ceramide_smiles)
    if ceramide is None:
        return False, "Invalid ceramide SMILES string"
    
    if not mol.HasSubstructMatch(ceramide):
        return False, "No ceramide group found"

    return True, "Molecule is a beta-D-glucosylceramide"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:83264',
                          'name': 'beta-D-glucosylceramide',
                          'definition': 'Any glucosylceramide in which the '
                                        'glucosyl component has '
                                        'beta-D-configuration.',
                          'parents': ['CHEBI:22798', 'CHEBI:36500']},
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
    'num_true_negatives': 14,
    'num_false_negatives': 14,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}