"""
Classifies: CHEBI:35436 D-glucoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_D_glucoside(smiles: str):
    """
    Determines if a molecule is a D-glucoside.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a D-glucoside, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the substructure pattern for D-glucose
    d_glucose_smiles = "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O"
    d_glucose_pattern = Chem.MolFromSmarts(d_glucose_smiles)
    
    if d_glucose_pattern is None:
        return False, "Unable to generate D-glucose pattern"

    # Check if the molecule contains the D-glucose substructure
    if not mol.HasSubstructMatch(d_glucose_pattern):
        return False, "Molecule does not contain D-glucose substructure"

    # Check for glycosidic bond (O-glycosidic linkage)
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O' and atom.GetDegree() == 2:
            neighbors = atom.GetNeighbors()
            if (neighbors[0].GetSymbol() == 'C' and neighbors[1].GetSymbol() == 'C' and
                mol.HasSubstructMatch(d_glucose_pattern, atom.GetIdx())):
                return True, "Molecule contains D-glucose substructure and glycosidic bond"

    return False, "Molecule does not contain glycosidic bond with D-glucose substructure"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35436',
                          'name': 'D-glucoside',
                          'definition': 'Any glucoside in which the glycoside '
                                        'group is derived from D-glucose.',
                          'parents': ['CHEBI:24278']},
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
    'num_true_positives': 79,
    'num_false_positives': 3,
    'num_true_negatives': 17,
    'num_false_negatives': 2,
    'precision': 0.9634146341463414,
    'recall': 0.9753086419753086,
    'f1': 0.9693251533742332,
    'accuracy': None}