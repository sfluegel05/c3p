"""
Classifies: CHEBI:35505 hydrate
"""
from rdkit import Chem

def is_hydrate(smiles: str):
    """
    Determines if a molecule is a hydrate (an addition compound that contains water in weak chemical combination with another compound).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydrate, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of water molecules
    water_smiles = "O"
    water_mol = Chem.MolFromSmiles(water_smiles)
    if water_mol is None:
        return False, "Error in water molecule definition"

    # Get the substructure matches for water in the molecule
    matches = mol.GetSubstructMatches(water_mol)
    if not matches:
        return False, "No water molecules found"

    # Check that there are additional atoms in the molecule besides the water molecules
    if mol.GetNumAtoms() == len(matches):
        return False, "Molecule only contains water"

    return True, "Molecule contains water in weak chemical combination with another compound"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35505',
                          'name': 'hydrate',
                          'definition': 'An addition compound that contains '
                                        'water in weak chemical combination '
                                        'with another compound.',
                          'parents': ['CHEBI:35504']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[23:31:44] WARNING: not removing hydrogen atom without '
             'neighbors\n'
             '[23:31:44] WARNING: not removing hydrogen atom without '
             'neighbors\n'
             '[23:31:44] WARNING: not removing hydrogen atom without '
             'neighbors\n',
    'stdout': '',
    'num_true_positives': 21,
    'num_false_positives': 19,
    'num_true_negatives': 1,
    'num_false_negatives': 0,
    'precision': 0.525,
    'recall': 1.0,
    'f1': 0.6885245901639345,
    'accuracy': None}