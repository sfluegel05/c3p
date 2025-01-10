"""
Classifies: CHEBI:33913 corrinoid
"""
"""
Classifies: CHEBI:33913 corrinoid
"""

from rdkit import Chem

def is_corrinoid(smiles: str):
    """
    Determines if a molecule is a corrinoid based on its SMILES string.
    A corrinoid is a derivative of the corrin nucleus, which contains four reduced or partly reduced
    pyrrole rings joined in a macrocycle by three =C- groups and one direct carbon-carbon bond linking
    alpha positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a corrinoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if molecule contains cobalt
    has_cobalt = any(atom.GetAtomicNum() == 27 for atom in mol.GetAtoms())
    if not has_cobalt:
        return False, "Molecule does not contain cobalt"

    # Find macrocycles (rings of size >= 13)
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    macrocycles = [ring for ring in atom_rings if len(ring) >= 13]

    if not macrocycles:
        return False, "No macrocycles of size >= 13 found"

    for ring in macrocycles:
        atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
        atomic_nums = [atom.GetAtomicNum() for atom in atoms_in_ring]

        # Count number of nitrogen and cobalt atoms in the ring
        num_N = atomic_nums.count(7)
        num_Co = atomic_nums.count(27)

        # Check if ring contains at least 4 nitrogen atoms and 1 cobalt atom
        if num_N >= 4 and num_Co >= 1:
            return True, "Molecule contains macrocycle with cobalt and nitrogen atoms characteristic of corrin nucleus"

    return False, "Molecule does not contain corrin nucleus"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:33913',
        'name': 'corrinoid',
        'definition': 'A derivative of the corrin nucleus, which contains four reduced or partly reduced pyrrole rings joined in a macrocycle by three =C- groups and one direct carbon-carbon bond linking alpha positions.',
        'parents': []
    },
    'config': {
        # Configuration parameters can be included here if necessary
    },
    'message': None,
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    # Metrics placeholders
    'num_true_positives': None,
    'num_false_positives': None,
    'num_true_negatives': None,
    'num_false_negatives': None,
    'num_negatives': None,
    'precision': None,
    'recall': None,
    'f1': None,
    'accuracy': None
}