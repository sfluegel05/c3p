"""
Classifies: CHEBI:26455 pyrroles
"""
from rdkit import Chem

def is_pyrroles(smiles: str):
    """
    Determines if a molecule is a pyrrole (An azole that includes only one N atom and no other heteroatom as a part of the aromatic skeleton).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pyrrole, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the aromatic ring information
    rings = mol.GetRingInfo()

    # Check for at least one 5-membered ring
    if not any(len(ring) == 5 for ring in rings.AtomRings()):
        return False, "No 5-membered rings found"

    # Find all aromatic 5-membered rings
    aromatic_rings = []
    for ring in rings.AtomRings():
        if len(ring) == 5:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if all(atom.GetIsAromatic() for atom in atoms):
                aromatic_rings.append(ring)

    if not aromatic_rings:
        return False, "No aromatic 5-membered rings found"

    # Check if the ring contains exactly one nitrogen atom and no other heteroatoms
    for ring in aromatic_rings:
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        nitrogen_count = sum(1 for atom in atoms if atom.GetSymbol() == 'N')
        other_heteroatom_count = sum(1 for atom in atoms if atom.GetSymbol() not in {'C', 'N', 'H'})

        if nitrogen_count == 1 and other_heteroatom_count == 0:
            return True, "Pyrrole structure found"

    return False, "No pyrrole structure found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26455',
                          'name': 'pyrroles',
                          'definition': 'An azole that includes only one N '
                                        'atom and no other heteroatom as a '
                                        'part of the aromatic skeleton.',
                          'parents': ['CHEBI:68452']},
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
    'num_true_positives': 29,
    'num_false_positives': 1,
    'num_true_negatives': 19,
    'num_false_negatives': 0,
    'precision': 0.9666666666666667,
    'recall': 1.0,
    'f1': 0.983050847457627,
    'accuracy': None}