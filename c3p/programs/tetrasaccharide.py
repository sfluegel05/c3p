"""
Classifies: CHEBI:50126 tetrasaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_tetrasaccharide(smiles: str):
    """
    Determines if a molecule is a tetrasaccharide (an oligosaccharide comprising four monomeric monosaccharide units).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrasaccharide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()

    # Check for at least four 5 or 6-membered rings (common size for monosaccharides)
    monosaccharide_rings = [ring for ring in rings.AtomRings() if len(ring) in [5, 6]]
    if len(monosaccharide_rings) < 4:
        return False, "Less than four monosaccharide rings found"

    # Check if all rings are connected to form a tetrasaccharide
    ring_atoms = set()
    for ring in monosaccharide_rings:
        ring_atoms.update(ring)

    bonds = mol.GetBonds()
    connected_rings = 0

    for bond in bonds:
        begin_atom = bond.GetBeginAtomIdx()
        end_atom = bond.GetEndAtomIdx()
        if begin_atom in ring_atoms and end_atom in ring_atoms:
            connected_rings += 1

    if connected_rings < 3:
        return False, "Rings are not connected to form a tetrasaccharide"

    return True, "Valid tetrasaccharide"

# Example usage:
# smiles = "OC[C@H]1O[C@H](O[C@H]2[C@@H](O)[C@H](O)[C@@H](CO)O[C@@H]2OC[C@H]2O[C@H](OC[C@H]3OC(O)[C@@H](O)[C@@H](O)[C@@H]3O)[C@@H](O)[C@@H](O)[C@@H]2O)[C@@H](O)[C@@H](O)[C@@H]1O"
# print(is_tetrasaccharide(smiles))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:50126',
                          'name': 'tetrasaccharide',
                          'definition': 'An oligosaccharide comprising four '
                                        'monomeric monosaccharide units.',
                          'parents': ['CHEBI:50699']},
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
    'num_true_positives': 6,
    'num_false_positives': 4,
    'num_true_negatives': 7,
    'num_false_negatives': 5,
    'precision': 0.6,
    'recall': 0.5454545454545454,
    'f1': 0.5714285714285713,
    'accuracy': None}