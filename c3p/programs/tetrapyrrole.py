"""
Classifies: CHEBI:26932 tetrapyrrole
"""
from rdkit import Chem

def is_tetrapyrrole(smiles: str):
    """
    Determines if a molecule is a tetrapyrrole.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrapyrrole, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()

    # Check for at least four 5-membered rings (pyrrole)
    if not any(len(ring) == 5 for ring in rings.AtomRings()):
        return False, "No 5-membered rings found"

    pyrrole_rings = []
    for ring in rings.AtomRings():
        if len(ring) == 5:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if any(atom.GetSymbol() == 'N' for atom in atoms):
                pyrrole_rings.append(ring)

    if len(pyrrole_rings) < 4:
        return False, "Less than four pyrrole rings found"

    # Check for the correct connectivity between pyrrole rings
    for i, ring in enumerate(pyrrole_rings):
        for j in range(i + 1, len(pyrrole_rings)):
            ring2 = pyrrole_rings[j]
            if not any(mol.GetBondBetweenAtoms(a1, a2) for a1 in ring for a2 in ring2):
                return False, "Pyrrole rings are not correctly connected"

    return True, "Molecule is a tetrapyrrole"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26932',
                          'name': 'tetrapyrrole',
                          'definition': 'A natural pigment containing four '
                                        'pyrrole rings joined by one-carbon '
                                        'units linking position 2 of one '
                                        'pyrrole ring to position 5 of the '
                                        'next.',
                          'parents': ['CHEBI:33833', 'CHEBI:38077']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 7-8: malformed \\N character escape (<string>, line 1)',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}