"""
Classifies: CHEBI:25477 naphthalenes
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_naphthalenes(smiles: str):
    """
    Determines if a molecule is a naphthalene (benzenoid aromatic compound having a skeleton composed of two ortho-fused benzene rings).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a naphthalene, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the aromatic ring information
    rings = mol.GetRingInfo()

    # Check for at least two 6-membered rings
    six_membered_rings = [ring for ring in rings.AtomRings() if len(ring) == 6]
    if len(six_membered_rings) < 2:
        return False, "Less than two 6-membered rings found"

    # Check for ortho-fused benzene rings
    for ring1 in six_membered_rings:
        for ring2 in six_membered_rings:
            if ring1 != ring2:
                shared_atoms = set(ring1).intersection(set(ring2))
                if len(shared_atoms) == 2:
                    return True, "Naphthalene structure found"

    return False, "No ortho-fused benzene rings found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:25477',
                          'name': 'naphthalenes',
                          'definition': 'Any benzenoid aromatic compound '
                                        'having a skeleton composed of two '
                                        'ortho-fused benzene rings.',
                          'parents': ['CHEBI:33836', 'CHEBI:36785']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 25-26: malformed \\N character escape (<string>, line '
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