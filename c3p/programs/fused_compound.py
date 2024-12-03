"""
Classifies: CHEBI:35293 fused compound
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_fused_compound(smiles: str):
    """
    Determines if a molecule is a fused compound (polycyclic compound with more than one ring sharing at least two adjacent atoms).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fused compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo().AtomRings()

    if len(rings) < 2:
        return False, "Less than two rings found"

    # Check for fused rings (rings sharing at least two adjacent atoms)
    for i, ring1 in enumerate(rings):
        for j, ring2 in enumerate(rings):
            if i >= j:
                continue
            common_atoms = set(ring1) & set(ring2)
            if len(common_atoms) >= 2:
                return True, "Fused rings found"

    return False, "No fused rings found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35293',
                          'name': 'fused compound',
                          'definition': 'A polycyclic compound that contains '
                                        'more than one ring with at least two '
                                        'common atoms (also known as '
                                        'bridgehead carbons) that are adjacent '
                                        'to each other.',
                          'parents': ['CHEBI:33635']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 28-29: malformed \\N character escape (<string>, line '
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