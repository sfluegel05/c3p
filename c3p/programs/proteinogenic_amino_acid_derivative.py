"""
Classifies: CHEBI:83811 proteinogenic amino acid derivative
"""
from rdkit import Chem

def is_proteinogenic_amino_acid_derivative(smiles: str):
    """
    Determines if a molecule is a proteinogenic amino acid derivative.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a proteinogenic amino acid derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS patterns for amino group, carboxy group, and side-chain functional groups
    amino_group = Chem.MolFromSmarts('[NX3;H2,H1;!$(NC=O)]')
    carboxy_group = Chem.MolFromSmarts('C(=O)[O;H1,H0-]')
    side_chain_functional_groups = [
        Chem.MolFromSmarts('[OH]'),
        Chem.MolFromSmarts('[SH]'),
        Chem.MolFromSmarts('[NH2]'),
        Chem.MolFromSmarts('[C=O]'),
        Chem.MolFromSmarts('[C#N]'),
        Chem.MolFromSmarts('[N+]'),
        Chem.MolFromSmarts('[S=O]'),
        Chem.MolFromSmarts('[P=O]')
    ]

    # Check for amino group
    if mol.HasSubstructMatch(amino_group):
        return True, "Contains amino group"

    # Check for carboxy group
    if mol.HasSubstructMatch(carboxy_group):
        return True, "Contains carboxy group"

    # Check for side-chain functional groups
    for group in side_chain_functional_groups:
        if mol.HasSubstructMatch(group):
            return True, "Contains side-chain functional group"

    # Check for replacement of any hydrogen by a heteroatom
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in [7, 8, 15, 16]:  # N, O, P, S
            if any(neighbor.GetAtomicNum() == 1 for neighbor in atom.GetNeighbors()):  # Check for hydrogen neighbors
                return True, "Contains hydrogen replaced by heteroatom"

    return False, "Does not match criteria for proteinogenic amino acid derivative"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:83811',
                          'name': 'proteinogenic amino acid derivative',
                          'definition': 'Any derivative of a proteinogenic '
                                        'amino acid resulting from reaction at '
                                        'an amino group, carboxy group, or a '
                                        'side-chain functional group, or from '
                                        'the replacement of any hydrogen by a '
                                        'heteroatom.',
                          'parents': ['CHEBI:83821']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': 'Python argument types in\n'
             '    Mol.HasSubstructMatch(Mol, NoneType)\n'
             'did not match C++ signature:\n'
             '    HasSubstructMatch(RDKit::ROMol self, RDKit::MolBundle query, '
             'RDKit::SubstructMatchParameters params=True)\n'
             '    HasSubstructMatch(RDKit::ROMol self, RDKit::ROMol query, '
             'RDKit::SubstructMatchParameters params)\n'
             '    HasSubstructMatch(RDKit::ROMol self, RDKit::MolBundle query, '
             'bool recursionPossible=True, bool useChirality=False, bool '
             'useQueryQueryMatches=False)\n'
             '    HasSubstructMatch(RDKit::ROMol self, RDKit::ROMol query, '
             'bool recursionPossible=True, bool useChirality=False, bool '
             'useQueryQueryMatches=False)',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}