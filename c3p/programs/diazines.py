"""
Classifies: CHEBI:38313 diazines
"""
from rdkit import Chem

def is_diazines(smiles: str):
    """
    Determines if a molecule is a diazine (any organic heterocyclic compound containing a benzene ring in which two of the C-H fragments have been replaced by isolobal nitrogens).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diazine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()

    # Check for at least one 6-membered ring
    if not any(len(ring) == 6 for ring in rings.AtomRings()):
        return False, "No 6-membered rings found"

    # Find all 6-membered rings
    six_membered_rings = [ring for ring in rings.AtomRings() if len(ring) == 6]

    for ring in six_membered_rings:
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        # Check if exactly two nitrogen atoms are present in the ring
        if sum(1 for atom in atoms if atom.GetSymbol() == 'N') == 2:
            # Ensure the remaining atoms are carbon
            if all(atom.GetSymbol() in ['C', 'N'] for atom in atoms):
                return True, "Molecule is a diazine"
    
    return False, "No 6-membered rings with exactly two nitrogen atoms found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:38313',
                          'name': 'diazines',
                          'definition': 'Any organic heterocyclic compound '
                                        'containing a benzene ring in which '
                                        'two of the C-H fragments have been '
                                        'replaced by isolobal nitrogens (the '
                                        'diazine parent structure).',
                          'parents': ['CHEBI:25693', 'CHEBI:38101']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 53-54: malformed \\N character escape (<string>, line '
             '1)',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}