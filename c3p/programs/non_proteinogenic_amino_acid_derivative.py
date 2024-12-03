"""
Classifies: CHEBI:83812 non-proteinogenic amino acid derivative
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_non_proteinogenic_amino_acid_derivative(smiles: str):
    """
    Determines if a molecule is a non-proteinogenic amino acid derivative.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a non-proteinogenic amino acid derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for amino group (NH2 or substituted)
    amino_group = False
    carboxy_group = False
    heteroatom_replacement = False

    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N':
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'H':
                    amino_group = True
                elif neighbor.GetSymbol() != 'C':
                    amino_group = True

        # Check for carboxy group (COOH or substituted)
        if atom.GetSymbol() == 'C':
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'O':
                    for next_neighbor in neighbor.GetNeighbors():
                        if next_neighbor.GetSymbol() == 'H':
                            carboxy_group = True
                        elif next_neighbor.GetSymbol() != 'C':
                            carboxy_group = True

        # Check for heteroatom replacement
        if atom.GetSymbol() not in ['C', 'H']:
            heteroatom_replacement = True

    if amino_group and carboxy_group and heteroatom_replacement:
        return True, "Molecule is a non-proteinogenic amino acid derivative"
    else:
        return False, "Molecule does not meet the criteria for non-proteinogenic amino acid derivative"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:83812',
                          'name': 'non-proteinogenic amino acid derivative',
                          'definition': 'Any derivative of a non-proteinogenic '
                                        'amino acid resulting from reaction at '
                                        'an amino group or carboxy group, or '
                                        'from the replacement of any hydrogen '
                                        'by a heteroatom.',
                          'parents': ['CHEBI:83821']},
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