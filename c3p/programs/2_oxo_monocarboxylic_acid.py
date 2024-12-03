"""
Classifies: CHEBI:35910 2-oxo monocarboxylic acid
"""
from rdkit import Chem

def is_2_oxo_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-oxo monocarboxylic acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (COOH)
    carboxylic_acid = False
    carboxylic_acid_carbon_idx = None
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetTotalDegree() == 3:
            neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
            if neighbors.count('O') == 2 and neighbors.count('C') == 1:
                carboxylic_acid = True
                carboxylic_acid_carbon_idx = atom.GetIdx()
                break
    
    if not carboxylic_acid:
        return False, "No carboxylic acid group found"

    # Check for 2-oxo group (C=O) adjacent to carboxylic acid
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetTotalDegree() == 3:
            neighbors = atom.GetNeighbors()
            if any(n.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2 for n in neighbors):
                for n in neighbors:
                    if n.GetSymbol() == 'C' and n.GetIdx() == carboxylic_acid_carbon_idx:
                        return True, "2-oxo monocarboxylic acid found"
    
    return False, "2-oxo group adjacent to carboxylic acid not found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35910',
                          'name': '2-oxo monocarboxylic acid',
                          'definition': 'Any monocarboxylic acid having a '
                                        '2-oxo substituent.',
                          'parents': ['CHEBI:35871']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 10,
    'num_false_positives': 1,
    'num_true_negatives': 10,
    'num_false_negatives': 1,
    'precision': 0.9090909090909091,
    'recall': 0.9090909090909091,
    'f1': 0.9090909090909091,
    'accuracy': None}