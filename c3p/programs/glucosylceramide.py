"""
Classifies: CHEBI:36500 glucosylceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_glucosylceramide(smiles: str):
    """
    Determines if a molecule is a glucosylceramide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glucosylceramide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the glucose moiety
    glucose_smiles = "C1(C(C(C(C(O1)CO)O)O)O)O"
    glucose_mol = Chem.MolFromSmiles(glucose_smiles)
    if not mol.HasSubstructMatch(glucose_mol):
        return False, "No glucose moiety found"

    # Check for the sphingosine backbone
    sphingosine_smiles = "C[C@@H](O)[C@H](NC=O)CO"
    sphingosine_mol = Chem.MolFromSmiles(sphingosine_smiles)
    if not mol.HasSubstructMatch(sphingosine_mol):
        return False, "No sphingosine backbone found"

    # Check for the ceramide structure
    ceramide_smiles = "C(CCCCCCCCCC)C(=O)N[C@@H](CO)C"
    ceramide_mol = Chem.MolFromSmiles(ceramide_smiles)
    if not mol.HasSubstructMatch(ceramide_mol):
        return False, "No ceramide structure found"

    return True, "Molecule is a glucosylceramide"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36500',
                          'name': 'glucosylceramide',
                          'definition': 'Any of the cerebrosides in which the '
                                        'monosaccharide head group is glucose.',
                          'parents': ['CHEBI:23079', 'CHEBI:62941']},
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
    'num_true_positives': 7,
    'num_false_positives': 8,
    'num_true_negatives': 6,
    'num_false_negatives': 7,
    'precision': 0.4666666666666667,
    'recall': 0.5,
    'f1': 0.4827586206896552,
    'accuracy': None}