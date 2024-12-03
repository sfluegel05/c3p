"""
Classifies: CHEBI:83820 non-proteinogenic amino acid
"""
from rdkit import Chem

def is_non_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a non-proteinogenic amino acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a non-proteinogenic amino acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for amino acid properties: presence of both amino group (NH2) and carboxyl group (COOH)
    has_amino_group = any(atom.GetSymbol() == 'N' and any(neigh.GetSymbol() == 'H' for neigh in atom.GetNeighbors()) for atom in mol.GetAtoms())
    has_carboxyl_group = any(atom.GetSymbol() == 'C' and any(neigh.GetSymbol() == 'O' for neigh in atom.GetNeighbors()) for atom in mol.GetAtoms())

    if not (has_amino_group and has_carboxyl_group):
        return False, "Molecule does not have both amino and carboxyl groups"

    # Check if it is a standard amino acid encoded in the genetic code
    standard_amino_acids_smiles = [
        'C(C(=O)O)N', 'CC(C(=O)O)N', 'CCC(C(=O)O)N', 'CC(C)C(C(=O)O)N', 'CC(C)CC(C(=O)O)N', 'CC(C(=O)O)N', 'C(C(=O)O)NCC(=O)O', 'C(C(=O)O)NCCC(=O)O', 'C(C(=O)O)NCCC(N)=O', 'C(C(=O)O)NCCS', 'C(C(=O)O)NCCO', 'C(C(=O)O)NCC(C)O', 'C(C(=O)O)NCC(C)C', 'C(C(=O)O)NCC(CO)O', 'C(C(=O)O)NCC(CO)O', 'C(C(=O)O)NCC(CO)O', 'C(C(=O)O)NCC(CO)O', 'C(C(=O)O)NCC(CO)O', 'C(C(=O)O)NCC(CO)O', 'C(C(=O)O)NCC(CO)O'
    ]

    if smiles in standard_amino_acids_smiles:
        return False, "Molecule is a standard amino acid encoded in the genetic code"

    return True, "Molecule is a non-proteinogenic amino acid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:83820',
                          'name': 'non-proteinogenic amino acid',
                          'definition': 'Any amino-acid that is not naturally '
                                        'encoded in the genetic code of any '
                                        'organism.',
                          'parents': ['CHEBI:33709']},
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
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 84,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}