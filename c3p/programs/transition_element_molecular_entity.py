"""
Classifies: CHEBI:33497 transition element molecular entity
"""
from rdkit import Chem


def is_transition_element_molecular_entity(smiles: str):
    """
    Determines if a molecule is a transition element molecular entity.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a transition element molecular entity, False otherwise
        str: Reason for classification
    """
    transition_elements = {
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
        "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
        "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn"
    }
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in transition_elements:
            return True, f"Contains transition element: {atom.GetSymbol()}"
    
    return False, "No transition elements found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33497',
                          'name': 'transition element molecular entity',
                          'definition': 'A molecular entity containing one or '
                                        'more atoms of a transition element.',
                          'parents': ['CHEBI:23367']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[23:02:45] Explicit valence for atom # 0 O, 3, is greater than '
             'permitted\n',
    'stdout': '',
    'num_true_positives': 74,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 2,
    'precision': 1.0,
    'recall': 0.9736842105263158,
    'f1': 0.9866666666666666,
    'accuracy': None}