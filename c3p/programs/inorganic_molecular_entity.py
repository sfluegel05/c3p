"""
Classifies: CHEBI:24835 inorganic molecular entity
"""
from rdkit import Chem

def is_inorganic_molecular_entity(smiles: str):
    """
    Determines if a molecule is an inorganic molecular entity (contains no carbon).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an inorganic molecular entity, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of carbon atoms
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            return False, "Contains carbon atoms"

    return True, "No carbon atoms found, classified as inorganic molecular entity"

# Examples
print(is_inorganic_molecular_entity("[Mn++].[O-]S([O-])(=O)=O"))  # manganese(II) sulfate
print(is_inorganic_molecular_entity("[O-2].[Cu+].[Cu+]"))  # cuprous oxide
print(is_inorganic_molecular_entity("[H][C]([H])([H])[H]"))  # methane (should return False)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:24835',
                          'name': 'inorganic molecular entity',
                          'definition': 'A molecular entity that contains no '
                                        'carbon.',
                          'parents': ['CHEBI:23367']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[20:41:52] Explicit valence for atom # 0 B, 5, is greater than '
             'permitted\n'
             '[20:41:52] Explicit valence for atom # 0 B, 5, is greater than '
             'permitted\n'
             '[20:41:52] WARNING: not removing hydrogen atom without '
             'neighbors\n'
             '[20:41:52] Explicit valence for atom # 0 B, 6, is greater than '
             'permitted\n'
             '[20:41:52] WARNING: not removing hydrogen atom without '
             'neighbors\n',
    'stdout': "(True, 'No carbon atoms found, classified as inorganic "
              "molecular entity')\n"
              "(True, 'No carbon atoms found, classified as inorganic "
              "molecular entity')\n"
              "(False, 'Contains carbon atoms')\n",
    'num_true_positives': 84,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 3,
    'precision': 1.0,
    'recall': 0.9655172413793104,
    'f1': 0.9824561403508771,
    'accuracy': None}