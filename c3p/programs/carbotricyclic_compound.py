"""
Classifies: CHEBI:38032 carbotricyclic compound
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_carbotricyclic_compound(smiles: str):
    """
    Determines if a molecule is a carbotricyclic compound (comprising three carbocyclic rings).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbotricyclic compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Check if there are at least three rings
    if len(atom_rings) < 3:
        return False, "Less than three rings found"

    # Check if all rings are carbocyclic (only contain carbon atoms)
    carbocyclic_rings = [ring for ring in atom_rings if all(mol.GetAtomWithIdx(atom_idx).GetSymbol() == 'C' for atom_idx in ring)]

    # Check if there are exactly three carbocyclic rings
    if len(carbocyclic_rings) != 3:
        return False, "Number of carbocyclic rings is not exactly three"

    return True, "Molecule is a carbotricyclic compound"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:38032',
                          'name': 'carbotricyclic compound',
                          'definition': 'A carbopolyclic compound comprising '
                                        'of three carbocyclic rings.',
                          'parents': ['CHEBI:35294', 'CHEBI:51959']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '[13:55:47] SMILES Parse Error: syntax error while parsing: '
             '[H]C(=O)C1=C/C[C@@]2([H])[C@]([H])(CC[C@]2(C)C[C@@]2([H])[C@]\x01([H])C(=O)C[C@@]2(C)O)[C@@H](C)CCC=C(C)C\n'
             '[13:55:47] SMILES Parse Error: Failed parsing SMILES '
             "'[H]C(=O)C1=C/C[C@@]2([H])[C@]([H])(CC[C@]2(C)C[C@@]2([H])[C@]\x01([H])C(=O)C[C@@]2(C)O)[C@@H](C)CCC=C(C)C' "
             'for input: '
             "'[H]C(=O)C1=C/C[C@@]2([H])[C@]([H])(CC[C@]2(C)C[C@@]2([H])[C@]\x01([H])C(=O)C[C@@]2(C)O)[C@@H](C)CCC=C(C)C'\n",
    'stdout': '',
    'num_true_positives': 98,
    'num_false_positives': 5,
    'num_true_negatives': 15,
    'num_false_negatives': 44,
    'precision': 0.9514563106796117,
    'recall': 0.6901408450704225,
    'f1': 0.8,
    'accuracy': None}