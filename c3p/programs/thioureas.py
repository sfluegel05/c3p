"""
Classifies: CHEBI:51276 thioureas
"""
from rdkit import Chem

def is_thioureas(smiles: str):
    """
    Determines if a molecule is a thiourea (compounds of general formula RR'NC(=S)NR''R''').

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiourea, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the thiourea substructure
    thiourea_smarts = 'NC(=S)N'
    thiourea_substructure = Chem.MolFromSmarts(thiourea_smarts)

    if mol.HasSubstructMatch(thiourea_substructure):
        return True, "Molecule contains the thiourea functional group"
    else:
        return False, "Molecule does not contain the thiourea functional group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:51276',
                          'name': 'thioureas',
                          'definition': 'Compounds of general formula '
                                        "RR'NC(=S)NR''R'''.",
                          'parents': ['CHEBI:50492']},
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
    'num_true_positives': 22,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 1,
    'precision': 1.0,
    'recall': 0.9565217391304348,
    'f1': 0.9777777777777777,
    'accuracy': None}