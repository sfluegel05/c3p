"""
Classifies: CHEBI:84144 L-phenylalanine derivative
"""
from rdkit import Chem

def is_L_phenylalanine_derivative(smiles: str):
    """
    Determines if a molecule is a L-phenylalanine derivative.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a L-phenylalanine derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMILES string for L-phenylalanine
    l_phenylalanine_smiles = "N[C@@H](Cc1ccccc1)C(=O)O"
    l_phenylalanine = Chem.MolFromSmiles(l_phenylalanine_smiles)

    # Check if the molecule contains the L-phenylalanine core structure
    if not mol.HasSubstructMatch(l_phenylalanine):
        return False, "Molecule does not contain the L-phenylalanine core structure"

    # Check for modifications at the amino group or carboxy group
    amino_group = Chem.MolFromSmarts("[N][C@@H](Cc1ccccc1)C(=O)O")
    carboxy_group = Chem.MolFromSmarts("N[C@@H](Cc1ccccc1)C(=O)[O]")

    if mol.HasSubstructMatch(amino_group) or mol.HasSubstructMatch(carboxy_group):
        return True, "Molecule is a L-phenylalanine derivative with modifications at amino or carboxy group"

    # Check for replacement of hydrogen by heteroatom
    heteroatoms = ["O", "N", "S", "P", "F", "Cl", "Br", "I"]
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in heteroatoms:
            return True, "Molecule is a L-phenylalanine derivative with heteroatom substitution"

    return False, "Molecule does not meet the criteria for L-phenylalanine derivative"

# Example usage:
# print(is_L_phenylalanine_derivative("COC1=CC=C(C[C@H](N)C(O)=O)C=C1O"))
# print(is_L_phenylalanine_derivative("N[C@@H](Cc1cc(I)c(O)c(I)c1)C(O)=O"))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:84144',
                          'name': 'L-phenylalanine derivative',
                          'definition': 'A proteinogenic amino acid derivative '
                                        'resulting from reaction of '
                                        'L-phenylalanine  at the amino group '
                                        'or the carboxy group, or from the '
                                        'replacement of any hydrogen of '
                                        'L-phenylalanine  by a heteroatom.',
                          'parents': ['CHEBI:25985', 'CHEBI:83811']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 10,
    'num_false_positives': 9,
    'num_true_negatives': 2,
    'num_false_negatives': 1,
    'precision': 0.5263157894736842,
    'recall': 0.9090909090909091,
    'f1': 0.6666666666666666,
    'accuracy': None}