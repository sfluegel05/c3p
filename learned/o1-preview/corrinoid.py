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

    # Define SMARTS pattern for corrin nucleus (approximate)
    # This pattern looks for four pyrrole-like rings connected in a macrocycle
    corrin_smarts = '[#7]-[#6]:[#6]-[#6]:[#7]-[#6]:[#6]-[#6]:[#7]-[#6]:[#6]-[#6]:[#7]-[#6]:[#6]-[#6]:[#6]-[#7]'

    corrin_pattern = Chem.MolFromSmarts(corrin_smarts)
    if corrin_pattern is None:
        return False, "Invalid SMARTS pattern for corrin nucleus"

    # Check if molecule has substructure match with corrin nucleus
    if mol.HasSubstructMatch(corrin_pattern):
        return True, "Molecule contains corrin nucleus"
    else:
        # As an alternative, check for ring systems with at least 4 nitrogen atoms
        ring_info = mol.GetRingInfo()
        atom_rings = ring_info.AtomRings()

        # Build ring systems by merging connected rings
        ring_systems = []
        for ring in atom_rings:
            ring_set = set(ring)
            merged = False
            for idx, existing_set in enumerate(ring_systems):
                if not ring_set.isdisjoint(existing_set):
                    ring_systems[idx] = existing_set.union(ring_set)
                    merged = True
                    break
            if not merged:
                ring_systems.append(ring_set)

        # Check each ring system
        for ring_sys in ring_systems:
            atoms_in_system = [mol.GetAtomWithIdx(idx) for idx in ring_sys]
            atomic_nums = [atom.GetAtomicNum() for atom in atoms_in_system]

            num_N = atomic_nums.count(7)
            num_atoms = len(atomic_nums)

            # Corrin nucleus typically has at least 15 atoms and 4 nitrogen atoms
            if num_N >= 4 and num_atoms >= 15:
                return True, "Molecule contains ring system characteristic of corrin nucleus"

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
    'attempt': 2,
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