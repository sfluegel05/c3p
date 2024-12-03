"""
Classifies: CHEBI:26822 sulfide
"""
from rdkit import Chem

def is_sulfide(smiles: str):
    """
    Determines if a molecule is a sulfide (Any sulfur molecular entity that involves either covalently bonded or anionic sulfur).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for sulfur atoms in the molecule
    sulfur_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'S']
    if not sulfur_atoms:
        return False, "No sulfur atoms found"

    # Check for covalent bonds or anionic sulfur
    for atom in sulfur_atoms:
        if atom.GetFormalCharge() == -1:
            return True, "Contains anionic sulfur"
        for neighbor in atom.GetNeighbors():
            if neighbor.GetSymbol() in ['C', 'H', 'O', 'N', 'S', 'P', 'F', 'Cl', 'Br', 'I']:
                return True, "Contains covalently bonded sulfur"

    return False, "No covalently bonded or anionic sulfur found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26822',
                          'name': 'sulfide',
                          'definition': 'Any sulfur molecular entity that '
                                        'involves either covalently bonded or '
                                        'anionic sulfur.',
                          'parents': ['CHEBI:26835']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[21:07:13] SMILES Parse Error: syntax error while parsing: '
             'COC1CC(=O)Nc2c(O)c(CC\\C=C(C)/C(O)C(C)C(C\\C=C\\C=C\\C=C\x01)OC(=O)C(C)NC(=O)C1CCCCC1)cc(O)c2SC\n'
             '[21:07:13] SMILES Parse Error: Failed parsing SMILES '
             "'COC1CC(=O)Nc2c(O)c(CC\\C=C(C)/C(O)C(C)C(C\\C=C\\C=C\\C=C\x01)OC(=O)C(C)NC(=O)C1CCCCC1)cc(O)c2SC' "
             'for input: '
             "'COC1CC(=O)Nc2c(O)c(CC\\C=C(C)/C(O)C(C)C(C\\C=C\\C=C\\C=C\x01)OC(=O)C(C)NC(=O)C1CCCCC1)cc(O)c2SC'\n"
             '[21:07:13] WARNING: not removing hydrogen atom without '
             'neighbors\n'
             '[21:07:13] WARNING: not removing hydrogen atom without '
             'neighbors\n',
    'stdout': '',
    'num_true_positives': 72,
    'num_false_positives': 3,
    'num_true_negatives': 17,
    'num_false_negatives': 2,
    'precision': 0.96,
    'recall': 0.972972972972973,
    'f1': 0.9664429530201343,
    'accuracy': None}