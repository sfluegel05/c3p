"""
Classifies: CHEBI:36785 carbobicyclic compound
"""
from rdkit import Chem

def is_carbobicyclic_compound(smiles: str):
    """
    Determines if a molecule is a carbobicyclic compound.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbobicyclic compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()

    # Check for exactly two rings
    if len(rings.AtomRings()) != 2:
        return False, "Molecule does not have exactly two rings"

    # Check if all ring atoms are carbon
    for ring in rings.AtomRings():
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() != 'C':
                return False, "Ring contains non-carbon atoms"

    # Check if the rings are fused (i.e., share at least one common atom)
    ring_atoms = [set(ring) for ring in rings.AtomRings()]
    if not ring_atoms[0].intersection(ring_atoms[1]):
        return False, "Rings are not fused"

    return True, "Molecule is a carbobicyclic compound"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36785',
                          'name': 'carbobicyclic compound',
                          'definition': 'A bicyclic compound in which all the '
                                        'ring atoms are carbon.',
                          'parents': ['CHEBI:33636', 'CHEBI:35294']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 1-2: malformed \\N character escape (<string>, line 1)',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}