"""
Classifies: CHEBI:36562 main-group coordination entity
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_main_group_coordination_entity(smiles: str):
    """
    Determines if a molecule is a main-group coordination entity.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a main-group coordination entity, False otherwise
        str: Reason for classification
    """
    main_group_elements = {'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'K', 'Ca', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Rb', 'Sr', 'In', 'Sn', 'Sb', 'Te', 'I', 'Cs', 'Ba', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn'}
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    atoms = mol.GetAtoms()
    central_atoms = [atom for atom in atoms if atom.GetSymbol() in main_group_elements and len(atom.GetNeighbors()) > 0]

    if not central_atoms:
        return False, "No central atom from main group elements found"

    for central_atom in central_atoms:
        neighbors = central_atom.GetNeighbors()
        if all(neighbor.GetSymbol() not in main_group_elements for neighbor in neighbors):
            return False, f"Central atom {central_atom.GetSymbol()} does not have ligands from main group elements"

    return True, "Main-group coordination entity found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36562',
                          'name': 'main-group coordination entity',
                          'definition': 'A coordination entity in which the '
                                        'central atom to which the ligands are '
                                        'attached comes from groups 1, 2, 13, '
                                        '14, 15, 16, 17, or 18 of the periodic '
                                        'table.',
                          'parents': ['CHEBI:33240', 'CHEBI:33579']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[23:55:56] Explicit valence for atom # 1 Te, 8, is greater than '
             'permitted\n'
             '[23:55:56] Explicit valence for atom # 1 Be, 6, is greater than '
             'permitted\n'
             '[23:55:57] Explicit valence for atom # 2 Si, 8, is greater than '
             'permitted\n'
             '[23:55:57] Explicit valence for atom # 0 C, 6, is greater than '
             'permitted\n',
    'stdout': '',
    'num_true_positives': 9,
    'num_false_positives': 9,
    'num_true_negatives': 4,
    'num_false_negatives': 4,
    'precision': 0.5,
    'recall': 0.6923076923076923,
    'f1': 0.5806451612903226,
    'accuracy': None}