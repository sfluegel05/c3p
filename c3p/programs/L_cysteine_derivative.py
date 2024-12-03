"""
Classifies: CHEBI:83824 L-cysteine derivative
"""
from rdkit import Chem

def is_L_cysteine_derivative(smiles: str):
    """
    Determines if a molecule is an L-cysteine derivative.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an L-cysteine derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMILES string of L-cysteine
    l_cysteine_smiles = "N[C@@H](CS)C(=O)O"
    l_cysteine_mol = Chem.MolFromSmiles(l_cysteine_smiles)
    if l_cysteine_mol is None:
        return None, None  # Should not happen, but added for safety

    # Check if the molecule contains the L-cysteine scaffold
    l_cysteine_substruct = mol.HasSubstructMatch(l_cysteine_mol)
    if not l_cysteine_substruct:
        return False, "Molecule does not contain the L-cysteine scaffold"

    # Identify the amino, carboxy, and thiol groups in the L-cysteine scaffold
    l_cysteine_atoms = mol.GetSubstructMatch(l_cysteine_mol)
    if not l_cysteine_atoms:
        return None, None  # Should not happen, but added for safety

    # Get the atom indices in the molecule corresponding to the L-cysteine scaffold
    amino_group_idx = l_cysteine_atoms[0]  # N
    thiol_group_idx = l_cysteine_atoms[1]  # S
    carboxy_group_idx = l_cysteine_atoms[3]  # C in carboxyl group

    # Check for modifications at these groups
    modified = False
    for atom_idx in [amino_group_idx, carboxy_group_idx, thiol_group_idx]:
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() not in l_cysteine_atoms:
                modified = True
                break

    if modified:
        return True, "Molecule is an L-cysteine derivative with modifications"
    else:
        return False, "Molecule is not modified from L-cysteine"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:83824',
                          'name': 'L-cysteine derivative',
                          'definition': 'A proteinogenic amino acid derivative '
                                        'resulting from the formal reaction of '
                                        'L-cysteine at the amino group, '
                                        'carboxy group, or thiol group, or '
                                        'from the replacement of any hydrogen '
                                        'of L-cysteine by a heteroatom.',
                          'parents': ['CHEBI:23509', 'CHEBI:83811']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 4,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 8,
    'num_false_positives': 1,
    'num_true_negatives': 10,
    'num_false_negatives': 3,
    'precision': 0.8888888888888888,
    'recall': 0.7272727272727273,
    'f1': 0.7999999999999999,
    'accuracy': None}