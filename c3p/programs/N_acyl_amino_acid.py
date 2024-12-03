"""
Classifies: CHEBI:51569 N-acyl-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_N_acyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acyl-amino acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acyl-amino acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxamide group (N-acyl group)
    carboxamide_found = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N':
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C':
                    for n_neighbor in neighbor.GetNeighbors():
                        if n_neighbor.GetSymbol() == 'O' and n_neighbor.GetIdx() != atom.GetIdx():
                            carboxamide_found = True
                            break
                if carboxamide_found:
                    break
        if carboxamide_found:
            break

    if not carboxamide_found:
        return False, "No N-acyl group found"

    # Check for amino acid backbone (alpha-amino acid)
    amino_acid_found = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetDegree() == 4:
            neighbors = atom.GetNeighbors()
            if any(neighbor.GetSymbol() == 'C' and neighbor.GetDegree() == 3 for neighbor in neighbors):
                if any(neighbor.GetSymbol() == 'N' for neighbor in neighbors):
                    if any(neighbor.GetSymbol() == 'O' for neighbor in neighbors):
                        amino_acid_found = True
                        break

    if not amino_acid_found:
        return False, "No amino acid backbone found"

    return True, "Molecule is an N-acyl-amino acid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:51569',
                          'name': 'N-acyl-amino acid',
                          'definition': 'A carboxamide resulting from the '
                                        'formal condensation of a carboxylic '
                                        'acid with the amino group of an amino '
                                        'acid.',
                          'parents': [   'CHEBI:33575',
                                         'CHEBI:37622',
                                         'CHEBI:83821']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 32-33: malformed \\N character escape (<string>, line '
             '1)',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}