"""
Classifies: CHEBI:38686 dihydroxyflavone
"""
from rdkit import Chem

def is_dihydroxyflavone(smiles: str):
    """
    Determines if a molecule is a dihydroxyflavone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dihydroxyflavone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule contains a flavone core structure
    flavone_smarts = '[O]=C1C=CC2=CC=CC=C2O1'
    flavone = Chem.MolFromSmarts(flavone_smarts)
    if not mol.HasSubstructMatch(flavone):
        return False, "Molecule does not contain a flavone core structure"

    # Check for hydroxy groups
    hydroxy_smarts = '[OH]'
    hydroxy = Chem.MolFromSmarts(hydroxy_smarts)
    hydroxy_matches = mol.GetSubstructMatches(hydroxy)

    if len(hydroxy_matches) < 2:
        return False, "Molecule does not have at least two hydroxy groups"

    # Check if the hydroxy groups are attached to the flavone core
    flavone_match = mol.GetSubstructMatch(flavone)
    flavone_atoms = set(flavone_match)
    hydroxy_on_flavone = 0

    for match in hydroxy_matches:
        if any(neighbor.GetIdx() in flavone_atoms for neighbor in mol.GetAtomWithIdx(match[0]).GetNeighbors()):
            hydroxy_on_flavone += 1

    if hydroxy_on_flavone < 2:
        return False, "Not all hydroxy groups are attached to the flavone core"

    return True, "Molecule is a dihydroxyflavone"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:38686',
                          'name': 'dihydroxyflavone',
                          'definition': 'Any hydroxyflavone in which two ring '
                                        'hydrogens are replaced by hydroxy '
                                        'substituents.',
                          'parents': ['CHEBI:24698']},
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
    'num_true_negatives': 11,
    'num_false_negatives': 11,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}