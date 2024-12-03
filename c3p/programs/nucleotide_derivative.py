"""
Classifies: CHEBI:231540 nucleotide derivative
"""
from rdkit import Chem

def is_nucleotide_derivative(smiles: str):
    """
    Determines if a molecule is a nucleotide derivative (a nucleoside phosphate that is derived from a nucleotide).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleotide derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a nucleoside phosphate group
    nucleoside_phosphate_smarts = Chem.MolFromSmarts("O[C@H]1[C@@H](O)[C@H](O[C@@H]1COP(=O)(O)OP(=O)(O)O)n1cnc2c1ncnc2N")
    if not mol.HasSubstructMatch(nucleoside_phosphate_smarts):
        nucleoside_phosphate_smarts = Chem.MolFromSmarts("O[C@H]1[C@@H](O)[C@H](O[C@@H]1COP(=O)(O)OP(=O)(O)O)n1cnc2c1ncnc2")
        if not mol.HasSubstructMatch(nucleoside_phosphate_smarts):
            return False, "No nucleoside phosphate group found"

    return True, "Molecule is a nucleotide derivative"

# Example usage
smiles = "CC(C)=CCC\\C(C)=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
result, reason = is_nucleotide_derivative(smiles)
print(result, reason)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:231540',
                          'name': 'nucleotide derivative',
                          'definition': 'A nucleoside phosphate that is '
                                        'derived from a nucleotide.',
                          'parents': ['CHEBI:25608']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 4,
    'success': True,
    'best': True,
    'error': '',
    'stdout': 'True Molecule is a nucleotide derivative\n',
    'num_true_positives': 82,
    'num_false_positives': 2,
    'num_true_negatives': 18,
    'num_false_negatives': 28,
    'precision': 0.9761904761904762,
    'recall': 0.7454545454545455,
    'f1': 0.845360824742268,
    'accuracy': None}