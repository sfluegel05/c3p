"""
Classifies: CHEBI:26400 purine ribonucleotide
"""
from rdkit import Chem

def is_purine_ribonucleotide(smiles: str):
    """
    Determines if a molecule is a purine ribonucleotide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a purine ribonucleotide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define purine base substructure
    purine_smarts = 'c1nc2[nH]cnc2c(=O)n1'
    purine_substructure = Chem.MolFromSmarts(purine_smarts)
    
    if not mol.HasSubstructMatch(purine_substructure):
        return False, "Molecule does not contain a purine base"

    # Define ribose substructure
    ribose_smarts = 'C1OC(O)C(O)C1O'
    ribose_substructure = Chem.MolFromSmarts(ribose_smarts)
    
    if not mol.HasSubstructMatch(ribose_substructure):
        return False, "Molecule does not contain a ribose sugar"

    # Check for phosphate group
    phosphate_smarts = 'P(=O)(O)O'
    phosphate_substructure = Chem.MolFromSmarts(phosphate_smarts)
    
    if not mol.HasSubstructMatch(phosphate_substructure):
        return False, "Molecule does not contain a phosphate group"

    return True, "Molecule is a purine ribonucleotide"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26400',
                          'name': 'purine ribonucleotide',
                          'definition': 'Any ribonucleotide that has a purine '
                                        'nucleobase.',
                          'parents': ['CHEBI:26395', 'CHEBI:26561']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 31-32: malformed \\N character escape (<string>, line '
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