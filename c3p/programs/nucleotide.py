"""
Classifies: CHEBI:36976 nucleotide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_nucleotide(smiles: str):
    """
    Determines if a molecule is a nucleotide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleotide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a nucleoside core (a sugar bonded to a nucleobase)
    nucleobases = ["c1ccn([C@H]2O[C@H](COP(O)(O)=O)[C@H](O)[C@H]2O)c(=O)n1",  # Cytosine
                   "n1cnc2c1ncnc2N",  # Adenine
                   "n1cnc2c1ncnc2O",  # Guanine
                   "n1ccc2c1ncnc2",  # Thymine
                   "n1ccncn1"]  # Uracil
    
    found_nucleobase = False
    for nucleobase in nucleobases:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(nucleobase)):
            found_nucleobase = True
            break

    if not found_nucleobase:
        return False, "No nucleobase found"

    # Check for the presence of a phosphate group
    phosphate = Chem.MolFromSmarts('P(=O)(O)(O)O')
    if not mol.HasSubstructMatch(phosphate):
        return False, "No phosphate group found"

    return True, "Molecule is a nucleotide"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36976',
                          'name': 'nucleotide',
                          'definition': 'A nucleotide is a nucleoside '
                                        'phosphate resulting from the '
                                        'condensation of the 3 or 5 hydroxy '
                                        'group of a nucleoside with phosphoric '
                                        'acid.',
                          'parents': ['CHEBI:25608']},
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
    'num_true_positives': 23,
    'num_false_positives': 19,
    'num_true_negatives': 1,
    'num_false_negatives': 20,
    'precision': 0.5476190476190477,
    'recall': 0.5348837209302325,
    'f1': 0.5411764705882354,
    'accuracy': None}