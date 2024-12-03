"""
Classifies: CHEBI:17855 triglyceride
"""
from rdkit import Chem

def is_triglyceride(smiles: str):
    """
    Determines if a molecule is a triglyceride (a glyceride resulting from the condensation of all three hydroxy groups of glycerol with fatty acids).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triglyceride, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone
    glycerol_backbone = Chem.MolFromSmarts("C(CO)CO")
    if not mol.HasSubstructMatch(glycerol_backbone):
        return False, "No glycerol backbone found"

    # Check for three ester linkages
    ester_linkage = Chem.MolFromSmarts("C(=O)O")
    ester_count = len(mol.GetSubstructMatches(ester_linkage))
    if ester_count != 3:
        return False, f"Expected 3 ester linkages, found {ester_count}"

    # Ensure each ester linkage is connected to the glycerol backbone
    glycerol_esters = Chem.MolFromSmarts("C(C(=O)O)C(C(=O)O)COC(=O)O")
    if not mol.HasSubstructMatch(glycerol_esters):
        return False, "Ester linkages are not connected to glycerol backbone properly"

    return True, "Molecule is a triglyceride"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:17855',
                          'name': 'triglyceride',
                          'definition': 'Any glyceride resulting from the '
                                        'condensation of all three hydroxy '
                                        'groups of glycerol '
                                        '(propane-1,2,3-triol) with fatty '
                                        'acids.',
                          'parents': ['CHEBI:47778', 'CHEBI:76579']},
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
    'num_true_negatives': 20,
    'num_false_negatives': 173,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}