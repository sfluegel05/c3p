"""
Classifies: CHEBI:25985 phenylalanine derivative
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phenylalanine_derivative(smiles: str):
    """
    Determines if a molecule is a phenylalanine derivative.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phenylalanine derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the phenylalanine core structure
    phenylalanine_smiles = "N[C@@H](Cc1ccccc1)C(=O)O"
    phenylalanine_mol = Chem.MolFromSmiles(phenylalanine_smiles)
    if phenylalanine_mol is None:
        return False, "Invalid phenylalanine core structure"

    # Check if the molecule contains the phenylalanine core
    if not mol.HasSubstructMatch(phenylalanine_mol):
        return False, "Molecule does not contain the phenylalanine core structure"

    # Check for modifications at the amino group, carboxy group, or replacement of hydrogen by heteroatom
    amino_group = Chem.MolFromSmarts("[N;R0]")
    carboxy_group = Chem.MolFromSmarts("[C](=O)[O;R0]")
    heteroatom = Chem.MolFromSmarts("[#6;R0]")

    if mol.HasSubstructMatch(amino_group) or mol.HasSubstructMatch(carboxy_group) or mol.HasSubstructMatch(heteroatom):
        return True, "Molecule is a phenylalanine derivative"
    else:
        return False, "No modifications at the amino group, carboxy group, or replacement of hydrogen by heteroatom"

# Example usage
smiles = "C([C@@H](C(OC)=O)N)C1=CC(O)=C(C=C1)O"  # melevodopa
result, reason = is_phenylalanine_derivative(smiles)
print(result, reason)  # Output: True, Molecule is a phenylalanine derivative


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:25985',
                          'name': 'phenylalanine derivative',
                          'definition': 'An amino acid derivative resulting '
                                        'from reaction of alanine at the amino '
                                        'group or the carboxy group, or from '
                                        'the replacement of any hydrogen of '
                                        'phenylalanine  by a heteroatom. The '
                                        'definition normally excludes peptides '
                                        'containing phenylalanine  residues.',
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
             'position 22-23: malformed \\N character escape (<string>, line '
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