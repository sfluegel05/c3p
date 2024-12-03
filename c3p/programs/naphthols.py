"""
Classifies: CHEBI:25392 naphthols
"""
from rdkit import Chem

def is_naphthols(smiles: str):
    """
    Determines if a molecule is a naphthol (hydroxynaphthalene derivative with a single hydroxy substituent).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a naphthol, False otherwise
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

    if len(aromatic_rings) < 2:
        return False, "Less than two aromatic 6-membered rings found"

    # Check for naphthalene structure (two fused aromatic rings)
    fused_rings = []
    for ring1 in aromatic_rings:
        for ring2 in aromatic_rings:
            if ring1 != ring2 and len(set(ring1).intersection(set(ring2))) >= 2:
                fused_rings.append((ring1, ring2))

    if not fused_rings:
        return False, "No fused aromatic rings found"

    # Check for exactly one hydroxy group
    hydroxy_count = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O' and len(atom.GetNeighbors()) == 1:
            neighbor = atom.GetNeighbors()[0]
            if neighbor.GetSymbol() == 'C' and neighbor.GetIsAromatic():
                hydroxy_count += 1

    if hydroxy_count != 1:
        return False, f"Number of hydroxy groups found: {hydroxy_count}"

    return True, "Valid naphthol"

# Example usage
smiles_list = [
    "O=C1C=C(CO)C=2[C@H]1[C@@]3(C=4C(=CC=C(C4C2)O)O[C@@H]3OC(C)C)C",
    "S1C(=N[C@@H](C1)[C@@H]2SC[C@@H]3N2C(=O)[C@@]4(N(C)[C@H]3SC4)C)C5=C(O)C=6C(=CC=CC6)C(=C5)OC",
    "O=C1OC2=C(C3=C(O)C4=C(OC)C=C(C=C4C=C3C[C@](C1)(O)C)OC)C(OC)=CC=5C2=C(O)C=6C(=O)C[C@@](O)(C)CC6C5",
    "O=C1C2=C(O)C3=C(OC)C(=CC=C3C(=C2[C@H](O)[C@H]4C1=CC[C@](O)(C)C4)OC)C",
    "O(C=1C2=C(O)C=CC(=C2C=CC1)C3=C(O)C4=C(OC)C=CC=C4C=C3)C",
    "Cc1cccc2c(cc(O)cc12)C(O)=O",
    "C12=C(C(=CC(=C1N)/N=N/C3=CC=C(C=C3)C4=CC=C(C=C4)/N=N/C5=CC(=C6C(=C5O)C=CC=C6)S(O)(=O)=O)S(O)(=O)=O)C=CC=C2",
    "COc1c(O)c(C(C)=O)c(COC(=O)C(\C)=C/C)c2[C@@H](C)[C@@H](O)C=Cc12",
    "S(SC=1C(=CC=C2C=CC=CC12)O)C=3C(=CC=C4C=CC=CC34)O",
    "O=C1C(=O)C=2C=3C(C(OC)=CC2O)=C(C)C=C(C3C1=O)O",
    "O=C1OC=2C(O)=CC(=C3C2C1=C(O)C=C3C)OCC=C(C)C",
    "OC=1C=2C(=C(C=CC2C)C)C=CC1C",
    "O1C=2C=3C(C=C(C2)C)=CC=C(C3C14O[C@@H](OC)C5=C4C=CC=C5OC)O",
    "C12=CC(=CC(=C1C(=C(C=C2)N)N=NC=3C=CC(=CC3S(O)(=O)=O)NC(=O)C)O)S(O)(=O)=O",
    "O=C1OC(C=2C=3C(C(O)=CC2CCC)=C(O)C=CC3)=CC(=C1)O",
    "Nc1ccc2c(O)cc(cc2c1)S(O)(=O)=O",
    "O=C1O[C@H](O)C=2C(O)=C(O)C(=C3C2C1=C(O)C=C3C)O",
    "O=C1C2=C(O)C(=CC=C2C(O)=C3C1=CC=4C(=O)[C@](O)(C)[C@@H]([C@H](C34)O)O)[C@@H]5O[C@@H]([C@@H](O)CC5)C",
    "O=C(CC=1C(O)=C(OC)C2=C(O)C=CC=C2C1)C"
]

for smiles in smiles_list:
    result, reason = is_naphthols(smiles)
    print(f"SMILES: {smiles}, Result: {result}, Reason: {reason}")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:25392',
                          'name': 'naphthols',
                          'definition': 'Any hydroxynaphthalene derivative '
                                        'that has a single hydroxy '
                                        'substituent.',
                          'parents': ['CHEBI:24727']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': 'SMILES: '
              'O=C1C=C(CO)C=2[C@H]1[C@@]3(C=4C(=CC=C(C4C2)O)O[C@@H]3OC(C)C)C, '
              'Result: False, Reason: Less than two aromatic 6-membered rings '
              'found\n'
              'SMILES: '
              'S1C(=N[C@@H](C1)[C@@H]2SC[C@@H]3N2C(=O)[C@@]4(N(C)[C@H]3SC4)C)C5=C(O)C=6C(=CC=CC6)C(=C5)OC, '
              'Result: True, Reason: Valid naphthol\n'
              'SMILES: '
              'O=C1OC2=C(C3=C(O)C4=C(OC)C=C(C=C4C=C3C[C@](C1)(O)C)OC)C(OC)=CC=5C2=C(O)C=6C(=O)C[C@@](O)(C)CC6C5, '
              'Result: False, Reason: Number of hydroxy groups found: 2\n'
              'SMILES: '
              'O=C1C2=C(O)C3=C(OC)C(=CC=C3C(=C2[C@H](O)[C@H]4C1=CC[C@](O)(C)C4)OC)C, '
              'Result: True, Reason: Valid naphthol\n'
              'SMILES: O(C=1C2=C(O)C=CC(=C2C=CC1)C3=C(O)C4=C(OC)C=CC=C4C=C3)C, '
              'Result: False, Reason: Number of hydroxy groups found: 2\n'
              'SMILES: Cc1cccc2c(cc(O)cc12)C(O)=O, Result: True, Reason: Valid '
              'naphthol\n'
              'SMILES: '
              'C12=C(C(=CC(=C1N)/N=N/C3=CC=C(C=C3)C4=CC=C(C=C4)/N=N/C5=CC(=C6C(=C5O)C=CC=C6)S(O)(=O)=O)S(O)(=O)=O)C=CC=C2, '
              'Result: True, Reason: Valid naphthol\n'
              'SMILES: '
              'COc1c(O)c(C(C)=O)c(COC(=O)C(\\C)=C/C)c2[C@@H](C)[C@@H](O)C=Cc12, '
              'Result: False, Reason: Less than two aromatic 6-membered rings '
              'found\n'
              'SMILES: S(SC=1C(=CC=C2C=CC=CC12)O)C=3C(=CC=C4C=CC=CC34)O, '
              'Result: False, Reason: Number of hydroxy groups found: 2\n'
              'SMILES: O=C1C(=O)C=2C=3C(C(OC)=CC2O)=C(C)C=C(C3C1=O)O, Result: '
              'False, Reason: Number of hydroxy groups found: 2\n'
              'SMILES: O=C1OC=2C(O)=CC(=C3C2C1=C(O)C=C3C)OCC=C(C)C, Result: '
              'False, Reason: Number of hydroxy groups found: 2\n'
              'SMILES: OC=1C=2C(=C(C=CC2C)C)C=CC1C, Result: True, Reason: '
              'Valid naphthol\n'
              'SMILES: '
              'O1C=2C=3C(C=C(C2)C)=CC=C(C3C14O[C@@H](OC)C5=C4C=CC=C5OC)O, '
              'Result: True, Reason: Valid naphthol\n'
              'SMILES: '
              'C12=CC(=CC(=C1C(=C(C=C2)N)N=NC=3C=CC(=CC3S(O)(=O)=O)NC(=O)C)O)S(O)(=O)=O, '
              'Result: True, Reason: Valid naphthol\n'
              'SMILES: O=C1OC(C=2C=3C(C(O)=CC2CCC)=C(O)C=CC3)=CC(=C1)O, '
              'Result: False, Reason: Number of hydroxy groups found: 4\n'
              'SMILES: Nc1ccc2c(O)cc(cc2c1)S(O)(=O)=O, Result: True, Reason: '
              'Valid naphthol\n'
              'SMILES: O=C1O[C@H](O)C=2C(O)=C(O)C(=C3C2C1=C(O)C=C3C)O, Result: '
              'False, Reason: Number of hydroxy groups found: 4\n'
              'SMILES: '
              'O=C1C2=C(O)C(=CC=C2C(O)=C3C1=CC=4C(=O)[C@](O)(C)[C@@H]([C@H](C34)O)O)[C@@H]5O[C@@H]([C@@H](O)CC5)C, '
              'Result: False, Reason: Less than two aromatic 6-membered rings '
              'found\n'
              'SMILES: O=C(CC=1C(O)=C(OC)C2=C(O)C=CC=C2C1)C, Result: False, '
              'Reason: Number of hydroxy groups found: 2\n',
    'num_true_positives': 8,
    'num_false_positives': 0,
    'num_true_negatives': 19,
    'num_false_negatives': 11,
    'precision': 1.0,
    'recall': 0.42105263157894735,
    'f1': 0.5925925925925926,
    'accuracy': None}