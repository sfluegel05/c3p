"""
Classifies: CHEBI:24315 glutamic acid derivative
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_glutamic_acid_derivative(smiles: str):
    """
    Determines if a molecule is a glutamic acid derivative.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glutamic acid derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the glutamic acid core structure
    glutamic_acid_smiles = "C(CC(=O)O)C(C(=O)O)N"
    glutamic_acid_mol = Chem.MolFromSmiles(glutamic_acid_smiles)
    if glutamic_acid_mol is None:
        return False, "Error in defining glutamic acid core structure"

    # Check if the molecule contains the glutamic acid core
    if not mol.HasSubstructMatch(glutamic_acid_mol):
        return False, "Molecule does not contain the glutamic acid core structure"

    # Check for modifications at the amino group or carboxy groups
    amino_group = Chem.MolFromSmarts("[NX3][CX4H2][CX4](C(=O)O)C(=O)O")
    carboxy_groups = Chem.MolFromSmarts("C(=O)O")

    amino_matches = mol.GetSubstructMatches(amino_group)
    carboxy_matches = mol.GetSubstructMatches(carboxy_groups)

    if not amino_matches and not carboxy_matches:
        return False, "No modifications at the amino group or carboxy groups"

    # Check for replacement of hydrogen by a heteroatom
    heteroatom_replacement = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in ['C', 'H']:
            heteroatom_replacement = True
            break

    if not heteroatom_replacement:
        return False, "No heteroatom replacement detected"

    return True, "Molecule is a glutamic acid derivative"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:24315',
                          'name': 'glutamic acid derivative',
                          'definition': 'An amino acid derivative resulting '
                                        'from reaction of glutamic acid at the '
                                        'amino group or either of the carboxy '
                                        'groups, or from the replacement of '
                                        'any hydrogen by a heteroatom. The '
                                        'definition normally excludes peptides '
                                        'containing glutamic acid residues.',
                          'parents': ['CHEBI:83821']},
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
    'num_true_positives': 12,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 12,
    'precision': 1.0,
    'recall': 0.5,
    'f1': 0.6666666666666666,
    'accuracy': None}