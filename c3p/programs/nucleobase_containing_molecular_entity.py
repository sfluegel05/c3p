"""
Classifies: CHEBI:61120 nucleobase-containing molecular entity
"""
from rdkit import Chem

def is_nucleobase_containing_molecular_entity(smiles: str):
    """
    Determines if a molecule is a nucleobase-containing molecular entity.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleobase-containing molecular entity, False otherwise
        str: Reason for classification
    """
    nucleobases = [
        'c1ncnc2n(cnc12)',  # Adenine
        'c1nc2ncnc2n1',  # Guanine
        'c1cnc[nH]c1',  # Cytosine
        'c1cc[nH]c(=O)n1',  # Thymine
        'c1cc[nH]c(=O)[nH]1',  # Uracil
    ]

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    for nucleobase in nucleobases:
        substructure = Chem.MolFromSmarts(nucleobase)
        if substructure is None:
            continue
        if mol.HasSubstructMatch(substructure):
            return True, f"Contains nucleobase: {nucleobase}"

    return False, "No nucleobase found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:61120',
                          'name': 'nucleobase-containing molecular entity',
                          'definition': 'Any compound that has a nucleobase as '
                                        'a part.',
                          'parents': ['CHEBI:33833', 'CHEBI:51143']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 3,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 195,
    'num_false_positives': 1,
    'num_true_negatives': 19,
    'num_false_negatives': 22,
    'precision': 0.9948979591836735,
    'recall': 0.8986175115207373,
    'f1': 0.9443099273607748,
    'accuracy': None}