"""
Classifies: CHEBI:33637 ortho-fused compound
"""
from rdkit import Chem

def is_ortho_fused_compound(smiles: str):
    """
    Determines if a molecule is an ortho-fused compound.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ortho-fused compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()
    atom_rings = rings.AtomRings()

    # Check for ortho-fused rings
    for i, ring1 in enumerate(atom_rings):
        for j, ring2 in enumerate(atom_rings):
            if i >= j:
                continue
            common_atoms = set(ring1).intersection(set(ring2))
            if len(common_atoms) == 2:
                return True, "Ortho-fused compound found"

    return False, "No ortho-fused rings found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33637',
                          'name': 'ortho-fused compound',
                          'definition': 'A polycyclic compound in which two '
                                        'rings have two, and only two, atoms '
                                        'in common. Such compounds have n '
                                        'common faces and 2n common atoms.',
                          'parents': ['CHEBI:35293']},
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
    'num_true_positives': 12,
    'num_false_positives': 7,
    'num_true_negatives': 5,
    'num_false_negatives': 0,
    'precision': 0.631578947368421,
    'recall': 1.0,
    'f1': 0.7741935483870968,
    'accuracy': None}