"""
Classifies: CHEBI:37948 oxaspiro compound
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_oxaspiro_compound(smiles: str):
    """
    Determines if a molecule is an oxaspiro compound (a spiro compound in which at least one of the cyclic components is an oxygen heterocycle).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oxaspiro compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    rings = mol.GetRingInfo()
    atom_rings = rings.AtomRings()
    bond_rings = rings.BondRings()

    spiro_atoms = set()
    for atom in mol.GetAtoms():
        if atom.GetDegree() == 4 and atom.IsInRing():
            spiro_atoms.add(atom.GetIdx())

    if not spiro_atoms:
        return False, "No spiro atoms found"

    for spiro_atom in spiro_atoms:
        connected_rings = [ring for ring in atom_rings if spiro_atom in ring]
        if len(connected_rings) < 2:
            continue

        for ring in connected_rings:
            ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if any(atom.GetSymbol() == 'O' for atom in ring_atoms):
                return True, "Found spiro atom with an oxygen-containing ring"

    return False, "No oxygen-containing ring found in spiro compound"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:37948',
                          'name': 'oxaspiro compound',
                          'definition': 'A spiro compound in which at least '
                                        'one of the cyclic components is an '
                                        'oxygen heterocyle.',
                          'parents': ['CHEBI:33599']},
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
    'num_true_positives': 27,
    'num_false_positives': 5,
    'num_true_negatives': 15,
    'num_false_negatives': 0,
    'precision': 0.84375,
    'recall': 1.0,
    'f1': 0.9152542372881356,
    'accuracy': None}