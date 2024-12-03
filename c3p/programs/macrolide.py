"""
Classifies: CHEBI:25106 macrolide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_macrolide(smiles: str):
    """
    Determines if a molecule is a macrolide (a macrocyclic lactone with a ring of twelve or more members derived from a polyketide).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a macrolide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all the rings in the molecule
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Check for macrocyclic lactone ring (12 or more members with at least one ester bond)
    for ring in atom_rings:
        if len(ring) >= 12:
            # Check if the ring contains an ester bond (C=O and O-C in the ring)
            contains_lactone = False
            for atom_idx in ring:
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetAtomicNum() == 8:  # Oxygen atom
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetIdx() in ring and neighbor.GetAtomicNum() == 6:  # Carbon atom
                            for carbon_neighbor in neighbor.GetNeighbors():
                                if carbon_neighbor.GetIdx() in ring and carbon_neighbor.GetAtomicNum() == 6:
                                    for bond in mol.GetBonds():
                                        if (bond.GetBeginAtomIdx() == neighbor.GetIdx() and bond.GetEndAtomIdx() == carbon_neighbor.GetIdx()) or \
                                           (bond.GetEndAtomIdx() == neighbor.GetIdx() and bond.GetBeginAtomIdx() == carbon_neighbor.GetIdx()):
                                            if bond.GetBondTypeAsDouble() == 2.0:  # Double bond (C=O)
                                                contains_lactone = True
                                                break
                            if contains_lactone:
                                break
                if contains_lactone:
                    break
            if contains_lactone:
                return True, "Macrocyclic lactone ring with 12 or more members found"

    return False, "No macrocyclic lactone ring with 12 or more members found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:25106',
                          'name': 'macrolide',
                          'definition': 'A macrocyclic lactone with a ring of '
                                        'twelve or more members derived from a '
                                        'polyketide.',
                          'parents': ['CHEBI:26188', 'CHEBI:63944']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[13:35:47] SMILES Parse Error: syntax error while parsing: '
             'CC[C@@H]1CC[C@@H]2O[C@@]3(O[C@@H](C[C@@H](C)O)[C@@H](C)CC3=O)[C@@H](C)[C@H](OC(=O)\\C=C\\[C@H](C)[C@@H](O)[C@H](C)C(=O)[C@H](C)[C@@H](O)[C@H](C)C(=O)[C@](C)(O)[C@@H](O)[C@H](C)C\\C=C\\C=C\x01)[C@@H]2C\n'
             '[13:35:47] SMILES Parse Error: Failed parsing SMILES '
             "'CC[C@@H]1CC[C@@H]2O[C@@]3(O[C@@H](C[C@@H](C)O)[C@@H](C)CC3=O)[C@@H](C)[C@H](OC(=O)\\C=C\\[C@H](C)[C@@H](O)[C@H](C)C(=O)[C@H](C)[C@@H](O)[C@H](C)C(=O)[C@](C)(O)[C@@H](O)[C@H](C)C\\C=C\\C=C\x01)[C@@H]2C' "
             'for input: '
             "'CC[C@@H]1CC[C@@H]2O[C@@]3(O[C@@H](C[C@@H](C)O)[C@@H](C)CC3=O)[C@@H](C)[C@H](OC(=O)\\C=C\\[C@H](C)[C@@H](O)[C@H](C)C(=O)[C@H](C)[C@@H](O)[C@H](C)C(=O)[C@](C)(O)[C@@H](O)[C@H](C)C\\C=C\\C=C\x01)[C@@H]2C'\n"
             '[13:35:47] SMILES Parse Error: unclosed ring for input: '
             "'[H][C@]12C=C(C)[C@@H](C)[C@@]3([H])[C@H](CC(C)C)NC(=O)[C@@]13OC(=O)\\C=C/[C@@H](O)[C@@H](O)[C@@H](C)C\\C(C)=C\x02'\n",
    'stdout': '',
    'num_true_positives': 2,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 154,
    'precision': 1.0,
    'recall': 0.01282051282051282,
    'f1': 0.02531645569620253,
    'accuracy': None}