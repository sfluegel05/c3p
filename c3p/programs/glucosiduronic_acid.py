"""
Classifies: CHEBI:24302 glucosiduronic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_glucosiduronic_acid(smiles: str):
    """
    Determines if a molecule is a glucosiduronic acid (glucuronic acid linked to another substance via a glycosidic bond).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glucosiduronic acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of glucuronic acid fragment
    glucuronic_acid_smiles = "O[C@@H]1[C@@H](O)[C@@H](O[C@@H]([C@H]1O)C(O)=O)"
    glucuronic_acid_mol = Chem.MolFromSmiles(glucuronic_acid_smiles)
    if glucuronic_acid_mol is None:
        return None, "Error in glucuronic acid SMILES string"

    # Check if the molecule contains the glucuronic acid fragment
    if not mol.HasSubstructMatch(glucuronic_acid_mol):
        return False, "No glucuronic acid fragment found"

    # Check for glycosidic bond (an oxygen atom linking glucuronic acid to another part of the molecule)
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O':
            neighbors = atom.GetNeighbors()
            if len(neighbors) == 2:
                if (neighbors[0].GetSymbol() == 'C' and neighbors[1].GetSymbol() == 'C') or (neighbors[0].GetSymbol() == 'C' and neighbors[1].GetSymbol() == 'H'):
                    return True, "Glucuronic acid linked via glycosidic bond"

    return False, "No glycosidic bond found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:24302',
                          'name': 'glucosiduronic acid',
                          'definition': 'Any substance produced by linking '
                                        'glucuronic acid to another substance '
                                        'via a glycosidic bond.',
                          'parents': ['CHEBI:35314', 'CHEBI:63436']},
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
    'num_false_negatives': 47,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}