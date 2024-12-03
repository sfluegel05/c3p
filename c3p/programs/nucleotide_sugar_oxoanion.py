"""
Classifies: CHEBI:59737 nucleotide-sugar oxoanion
"""
from rdkit import Chem

def is_nucleotide_sugar_oxoanion(smiles: str):
    """
    Determines if a molecule is a nucleotide-sugar oxoanion.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleotide-sugar oxoanion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for nucleotide-sugar backbone
    nucleotides = ["A", "C", "G", "T", "U"]
    nucleotide_found = any(nucleotide in smiles for nucleotide in nucleotides)
    if not nucleotide_found:
        return False, "No nucleotide found in the molecule"

    # Check for sugar moiety (common sugar identifiers in SMILES)
    sugars = ["C1OC", "C1CC", "C1C[C@", "C1C[C@@", "C1O[C@", "C1O[C@@", "C1CO"]
    sugar_found = any(sugar in smiles for sugar in sugars)
    if not sugar_found:
        return False, "No sugar moiety found in the molecule"

    # Check for deprotonated diphosphate groups (presence of [O-])
    deprotonated_phosphate = "[O-]"
    if deprotonated_phosphate not in smiles:
        return False, "No deprotonated diphosphate group found in the molecule"

    return True, "Molecule is a nucleotide-sugar oxoanion"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:59737',
                          'name': 'nucleotide-sugar oxoanion',
                          'definition': 'Any nucleotide-sugar in which one or '
                                        'more of the diphosphate OH groups has '
                                        'been deprotonated.',
                          'parents': ['CHEBI:58945']},
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
    'num_false_negatives': 26,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}