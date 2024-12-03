"""
Classifies: CHEBI:61293 adenyl nucleotide
"""
from rdkit import Chem

def is_adenyl_nucleotide(smiles: str):
    """
    Determines if a molecule is an adenyl nucleotide (a nucleotide having adenine as the base).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an adenyl nucleotide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for adenine base
    adenine_smarts = "n1cnc2c1ncnc2N"
    adenine = Chem.MolFromSmarts(adenine_smarts)
    if adenine is None:
        return False, "Error creating adenine molecule"

    # Check if adenine substructure is present in the molecule
    if not mol.HasSubstructMatch(adenine):
        return False, "Adenine base not found"

    # Check for the presence of ribose (sugar) and phosphate groups
    ribose_smarts = "C1(C(C(C(O1)CO)O)O)"
    phosphate_smarts = "P(=O)(O)O"

    ribose = Chem.MolFromSmarts(ribose_smarts)
    phosphate = Chem.MolFromSmarts(phosphate_smarts)

    if ribose is None or phosphate is None:
        return False, "Error creating ribose or phosphate molecule"

    if not mol.HasSubstructMatch(ribose):
        return False, "Ribose (sugar) not found"

    if not mol.HasSubstructMatch(phosphate):
        return False, "Phosphate group not found"

    return True, "Adenyl nucleotide"

# Example usage
smiles = "Nc1ncnc2n(cnc12)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OP(O)(O)=O)[C@@H](O)[C@H]1O"
print(is_adenyl_nucleotide(smiles))  # Should return True, "Adenyl nucleotide"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:61293',
                          'name': 'adenyl nucleotide',
                          'definition': 'A  nucleotide having adenine as the '
                                        'base.',
                          'parents': ['CHEBI:26395']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '',
    'stdout': "(True, 'Adenyl nucleotide')\n",
    'num_true_positives': 10,
    'num_false_positives': 3,
    'num_true_negatives': 7,
    'num_false_negatives': 0,
    'precision': 0.7692307692307693,
    'recall': 1.0,
    'f1': 0.8695652173913044,
    'accuracy': None}