"""
Classifies: CHEBI:35496 fluorobenzenes
"""
from rdkit import Chem

def is_fluorobenzenes(smiles: str):
    """
    Determines if a molecule is a fluorobenzene (benzene or substituted benzene carrying at least one fluoro group).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fluorobenzene, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the aromatic ring information
    rings = mol.GetRingInfo()

    # Check for at least one 6-membered ring
    if not any(len(ring) == 6 for ring in rings.AtomRings()):
        return False, "No 6-membered rings found"

    # Find all aromatic 6-membered rings
    aromatic_rings = []
    for ring in rings.AtomRings():
        if len(ring) == 6:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if all(atom.GetIsAromatic() for atom in atoms):
                aromatic_rings.append(ring)

    if not aromatic_rings:
        return False, "No aromatic 6-membered rings found"

    # Check if at least one aromatic ring contains a fluorine atom
    for ring in aromatic_rings:
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if any(atom.GetSymbol() == 'F' for atom in atoms):
            return True, "Molecule is a fluorobenzene"

    # Check substituents
    ring_atoms = set(aromatic_rings[0])
    for atom_idx in ring_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetSymbol() == 'F':
                return True, "Molecule is a fluorobenzene"

    return False, "No fluoro group found in the benzene ring"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35496',
                          'name': 'fluorobenzenes',
                          'definition': 'Any fluoroarene that is a benzene or '
                                        'a substituted benzene carrying at '
                                        'least one fluoro group.',
                          'parents': ['CHEBI:22712', 'CHEBI:37143']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 15-16: malformed \\N character escape (<string>, line '
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