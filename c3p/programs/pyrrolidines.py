"""
Classifies: CHEBI:38260 pyrrolidines
"""
from rdkit import Chem

def is_pyrrolidines(smiles: str):
    """
    Determines if a molecule is a pyrrolidine (heterocyclic amines having a saturated five-membered ring).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pyrrolidine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    ring_info = mol.GetRingInfo()

    # Check for at least one 5-membered ring
    if not any(len(ring) == 5 for ring in ring_info.AtomRings()):
        return False, "No 5-membered rings found"

    # Find all 5-membered rings
    five_membered_rings = [ring for ring in ring_info.AtomRings() if len(ring) == 5]

    for ring in five_membered_rings:
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        # Check if the ring is saturated (no double or triple bonds)
        if all(atom.GetDegree() == 2 or atom.GetDegree() == 3 for atom in atoms):
            # Check if the ring contains at least one nitrogen atom
            if any(atom.GetSymbol() == 'N' for atom in atoms):
                return True, "Saturated five-membered ring with nitrogen found"
    
    return False, "No saturated five-membered ring with nitrogen found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:38260',
                          'name': 'pyrrolidines',
                          'definition': 'Any of a class of heterocyclic amines '
                                        'having a saturated five-membered '
                                        'ring.',
                          'parents': ['CHEBI:25693', 'CHEBI:38101']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[00:14:00] SMILES Parse Error: syntax error while parsing: '
             'O=C1[C@@H](O)C[C@H]2[C@H](/C(/O)=C\x03/C(=O)N[C@H](C3=O)C)[C@H](C)C=C[C@@H]2C1\n'
             '[00:14:00] SMILES Parse Error: Failed parsing SMILES '
             "'O=C1[C@@H](O)C[C@H]2[C@H](/C(/O)=C\x03/C(=O)N[C@H](C3=O)C)[C@H](C)C=C[C@@H]2C1' "
             'for input: '
             "'O=C1[C@@H](O)C[C@H]2[C@H](/C(/O)=C\x03/C(=O)N[C@H](C3=O)C)[C@H](C)C=C[C@@H]2C1'\n",
    'stdout': '',
    'num_true_positives': 65,
    'num_false_positives': 2,
    'num_true_negatives': 18,
    'num_false_negatives': 3,
    'precision': 0.9701492537313433,
    'recall': 0.9558823529411765,
    'f1': 0.962962962962963,
    'accuracy': None}