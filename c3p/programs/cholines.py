"""
Classifies: CHEBI:23217 cholines
"""
from rdkit import Chem

def is_cholines(smiles: str):
    """
    Determines if a molecule is a choline or its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a choline or its derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for quaternary ammonium ion with choline structure
    quaternary_ammonium = False
    phosphate_group = False

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1:  # Nitrogen with +1 charge
            neighbors = [n.GetAtomicNum() for n in atom.GetNeighbors()]
            if neighbors.count(6) == 3 and neighbors.count(1) == 1:  # 3 carbons and 1 hydrogen
                quaternary_ammonium = True
        if atom.GetAtomicNum() == 15 and atom.GetFormalCharge() == 0:  # Phosphorus with no charge
            neighbors = [n.GetAtomicNum() for n in atom.GetNeighbors()]
            if neighbors.count(8) == 4:  # 4 oxygens
                phosphate_group = True

    if quaternary_ammonium and phosphate_group:
        return True, "Molecule is a choline derivative"
    else:
        return False, "Molecule does not match choline derivative structure"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:23217',
                          'name': 'cholines',
                          'definition': 'A quaternary ammonium ion based on '
                                        'the choline ion and its substituted '
                                        'derivatives thereof.',
                          'parents': ['CHEBI:35267']},
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
    'num_false_negatives': 216,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}