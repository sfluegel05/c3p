"""
Classifies: CHEBI:27116 trihydroxyflavone
"""
from rdkit import Chem

def is_trihydroxyflavone(smiles: str):
    """
    Determines if a molecule is a trihydroxyflavone (hydroxyflavone carrying three hydroxy groups at unspecified positions).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trihydroxyflavone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the flavone core structure (2-phenylchromen-4-one)
    flavone_pattern = Chem.MolFromSmarts('c1cc2c(cc1)oc(=O)c3c2cc(oc3)c4ccccc4')
    if not mol.HasSubstructMatch(flavone_pattern):
        return False, "Molecule does not contain the flavone core structure"

    # Count the number of hydroxy groups (-OH)
    hydroxy_count = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O':
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'H':
                    hydroxy_count += 1

    if hydroxy_count != 3:
        return False, f"Molecule contains {hydroxy_count} hydroxy groups, expected 3"

    return True, "Molecule is a trihydroxyflavone"

# Example usage:
# smiles = "CC(C)=CCc1c(oc2cc3OC(C)(C)C=Cc3c(O)c2c1=O)-c1ccc(O)cc1O"
# print(is_trihydroxyflavone(smiles))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:27116',
                          'name': 'trihydroxyflavone',
                          'definition': 'Any  hydroxyflavone carrying three '
                                        'hydroxy groups at unspecified '
                                        'positions.',
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
    'num_true_negatives': 15,
    'num_false_negatives': 15,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}