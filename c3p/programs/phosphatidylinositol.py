"""
Classifies: CHEBI:28874 phosphatidylinositol
"""
from rdkit import Chem

def is_phosphatidylinositol(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylinositol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of the inositol ring
    inositol_smiles = "C1(C(C(C(C(C1O)O)O)O)O)O"
    inositol_mol = Chem.MolFromSmiles(inositol_smiles)
    if not mol.HasSubstructMatch(inositol_mol):
        return False, "No inositol ring found"

    # Check for the presence of the glycerophospho group
    glycerophospho_smiles = "OCC(O)COP(O)(=O)O"
    glycerophospho_mol = Chem.MolFromSmiles(glycerophospho_smiles)
    if not mol.HasSubstructMatch(glycerophospho_mol):
        return False, "No glycerophospho group found"

    # Check for esterified phosphatidyl group
    ester_smiles = "C(=O)O"
    ester_mol = Chem.MolFromSmiles(ester_smiles)
    if not mol.HasSubstructMatch(ester_mol):
        return False, "No esterified phosphatidyl group found"

    return True, "Molecule is a phosphatidylinositol"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:28874',
                          'name': 'phosphatidylinositol',
                          'definition': 'Any glycerophosphoinositol having one '
                                        'phosphatidyl group esterified to one '
                                        'of the hydroxy groups of inositol.',
                          'parents': ['CHEBI:36315']},
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
    'num_true_positives': 48,
    'num_false_positives': 3,
    'num_true_negatives': 17,
    'num_false_negatives': 0,
    'precision': 0.9411764705882353,
    'recall': 1.0,
    'f1': 0.9696969696969697,
    'accuracy': None}