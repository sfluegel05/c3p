"""
Classifies: CHEBI:46874 alpha-amino acid ester
"""
from rdkit import Chem

def is_alpha_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid ester.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino acid ester, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of an alpha-amino acid backbone
    alpha_amino_acid = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetDegree() == 4:
            neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
            if 'N' in neighbors and 'C' in neighbors and 'C' in neighbors:
                alpha_amino_acid = True
                break

    if not alpha_amino_acid:
        return False, "No alpha-amino acid backbone found"

    # Check for ester group
    ester_group = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetDegree() == 3:
            neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
            if neighbors.count('O') == 2 and 'C' in neighbors:
                ester_group = True
                break

    if not ester_group:
        return False, "No ester group found"

    return True, "Alpha-amino acid ester found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:46874',
                          'name': 'alpha-amino acid ester',
                          'definition': 'The amino acid ester derivative '
                                        'obtained the formal condensation of '
                                        'an alpha-amino acid with an alcohol.',
                          'parents': ['CHEBI:46668']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 7,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 21,
    'precision': 1.0,
    'recall': 0.25,
    'f1': 0.4,
    'accuracy': None}