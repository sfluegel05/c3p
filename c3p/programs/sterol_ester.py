"""
Classifies: CHEBI:35915 sterol ester
"""
"""
Classifies: sterol ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.
    A sterol ester is a steroid ester obtained by formal condensation of the carboxy group of any carboxylic acid with the 3-hydroxy group of a sterol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Find rings of size 5 and 6
    ring_sizes = [len(r) for r in atom_rings]
    six_membered_rings = [r for r in atom_rings if len(r) == 6]
    five_membered_rings = [r for r in atom_rings if len(r) == 5]

    if len(six_membered_rings) < 3 or len(five_membered_rings) < 1:
        return False, "Does not contain required number of rings for sterol nucleus"

    # Build ring adjacency map
    # The ring adjacency can be determined by checking the bonds shared between rings
    ring_bond_map = {}
    for idx1, ring1 in enumerate(atom_rings):
        bonds1 = set()
        for i in range(len(ring1)):
            a1 = ring1[i]
            a2 = ring1[(i+1)%len(ring1)]
            bond = mol.GetBondBetweenAtoms(a1, a2).GetIdx()
            bonds1.add(bond)
        for idx2, ring2 in enumerate(atom_rings):
            if idx1 >= idx2:
                continue
            bonds2 = set()
            for i in range(len(ring2)):
                a1 = ring2[i]
                a2 = ring2[(i+1)%len(ring2)]
                bond = mol.GetBondBetweenAtoms(a1, a2).GetIdx()
                bonds2.add(bond)
            shared_bonds = bonds1.intersection(bonds2)
            if shared_bonds:
                ring_bond_map.setdefault(idx1, set()).add(idx2)
                ring_bond_map.setdefault(idx2, set()).add(idx1)

    # Identify fused ring system of size 6-6-6-5
    # Try to find a set of four rings that are connected and have sizes 6-6-6-5
    found_steroid_nucleus = False
    for idx1 in range(len(atom_rings)):
        for idx2 in ring_bond_map.get(idx1, []):
            for idx3 in ring_bond_map.get(idx2, []):
                for idx4 in ring_bond_map.get(idx3, []):
                    ring_sizes_sequence = [len(atom_rings[idx]) for idx in [idx1, idx2, idx3, idx4]]
                    if sorted(ring_sizes_sequence) == [5,6,6,6]:
                        unique_indices = set([idx1, idx2, idx3, idx4])
                        if len(unique_indices) == 4:
                            # Check if they are all connected
                            connections = [
                                idx2 in ring_bond_map[idx1],
                                idx3 in ring_bond_map[idx2],
                                idx4 in ring_bond_map[idx3],
                            ]
                            if all(connections):
                                found_steroid_nucleus = True
                                steroid_ring_indices = [idx1, idx2, idx3, idx4]
                                break
                if found_steroid_nucleus:
                    break
            if found_steroid_nucleus:
                break
        if found_steroid_nucleus:
            break

    if not found_steroid_nucleus:
        return False, "No steroid nucleus found"

    # Get atoms in steroid nucleus
    steroid_atom_indices = set()
    for idx in steroid_ring_indices:
        steroid_atom_indices.update(atom_rings[idx])

    # Find oxygen atoms attached to steroid nucleus
    steroid_atoms = [mol.GetAtomWithIdx(idx) for idx in steroid_atom_indices]
    steroid_atom_indices_set = set(steroid_atom_indices)
    oxygen_atoms = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() in steroid_atom_indices_set:
                    oxygen_atoms.append(atom)
                    break

    if not oxygen_atoms:
        return False, "No oxygen atoms attached to steroid nucleus"

    # Check if any of these oxygen atoms are part of an ester group
    ester_pattern = Chem.MolFromSmarts('[CX3](=O)[OX2H0]')  # Ester group
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    ester_oxygen_indices = [match[2] for match in ester_matches]  # index of oxygen atom in ester

    for o_atom in oxygen_atoms:
        if o_atom.GetIdx() in ester_oxygen_indices:
            return True, "Contains steroid nucleus with esterified oxygen attached"

    return False, "No esterified oxygen attached to steroid nucleus found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35915',
                          'name': 'sterol ester',
                          'definition': 'A steroid ester obtained by formal '
                                        'condensation of the carboxy group of '
                                        'any carboxylic acid with the '
                                        '3-hydroxy group of a sterol.',
                          'parents': ['CHEBI:33308', 'CHEBI:47880'],
                          'xrefs': ['KEGG:C01958'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'CNC(=O)CC[C@H](N)C(O)=O',
                                     'name': 'N(5)-methyl-L-glutamine',
                                     'reason': 'Does not contain required '
                                               'number of rings for sterol '
                                               'nucleus'},
                                 {   'smiles': 'C[C@H](OP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c1nc(N)[nH]c2=O)C(O)=O',
                                     'name': "L-lactyl-2-diphospho-5'-guanosine",
                                     'reason': 'Does not contain required '
                                               'number of rings for sterol '
                                               'nucleus'},
                                 {   'smiles': '[H][C@@]1(CC[C@]2(C)[C@]1([H])CC[C@]1([H])[C@@]3(C)CCCC(C)(C)[C@]3([H])CC[C@@]21C)C(=C)CCC=C(C)C',
                                     'name': 'dammara-20,24-diene',
                                     'reason': 'No oxygen atoms attached to '
                                               'steroid nucleus'},
                                 {   'smiles': 'O=C1C2=C(O)C3=C(O)C=CC=C3C(=C2C(C(=O)C)C(C1)(O)C)C=4C=5C(C(O)=C6C4C(C(=O)C)C(O)(C)CC6=O)=C(O)C=CC5',
                                     'name': 'A 39183A',
                                     'reason': 'Does not contain required '
                                               'number of rings for sterol '
                                               'nucleus'},
                                 {   'smiles': 'C[C@H](O)[C@@H](O)C1=Nc2c(NC1)[nH]c(N)nc2=O',
                                     'name': 'L-threo-7,8-dihydrobiopterin',
                                     'reason': 'Does not contain required '
                                               'number of rings for sterol '
                                               'nucleus'},
                                 {   'smiles': 'C[C@@H]1CN([C@H](COC2=C(C=C(C=C2)NC(=O)C)C(=O)N(C[C@H]1OC)C)C)CC3=CC=C(C=C3)F',
                                     'name': 'N-[(4S,7R,8S)-5-[(4-fluorophenyl)methyl]-8-methoxy-4,7,10-trimethyl-11-oxo-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]acetamide',
                                     'reason': 'Does not contain required '
                                               'number of rings for sterol '
                                               'nucleus'},
                                 {   'smiles': 'C(CCN1CCCCC1)(CCN(C(C)C)C(C)=O)(C(N)=O)C2=C(Cl)C=CC=C2',
                                     'name': 'bidisomide',
                                     'reason': 'Does not contain required '
                                               'number of rings for sterol '
                                               'nucleus'},
                                 {   'smiles': 'OCCCCCCCCCCC#CC=C',
                                     'name': '13-Tetradece-11-yn-1-ol',
                                     'reason': 'Does not contain required '
                                               'number of rings for sterol '
                                               'nucleus'},
                                 {   'smiles': 'C[C@@H]1CCCCO[C@@H]([C@@H](CN(C(=O)C2=C(O1)C=CC(=C2)N(C)C)[C@@H](C)CO)C)CN(C)CC3=CC=C(C=C3)OC',
                                     'name': '(3R,9S,10R)-16-(dimethylamino)-12-[(2S)-1-hydroxypropan-2-yl]-9-[[(4-methoxyphenyl)methyl-methylamino]methyl]-3,10-dimethyl-2,8-dioxa-12-azabicyclo[12.4.0]octadeca-1(14),15,17-trien-13-one',
                                     'reason': 'Does not contain required '
                                               'number of rings for sterol '
                                               'nucleus'},
                                 {   'smiles': 'O(CCC12CC3CC(C1)CC(C2)C3)C(=O)CC4=CC=C(OCC(O)CNC(C)C)C=C4',
                                     'name': 'Adaprolol',
                                     'reason': 'Does not contain required '
                                               'number of rings for sterol '
                                               'nucleus'}],
    'sample_false_negatives': [   {   'smiles': '[H][C@]12CC[C@]3(C)[C@]([H])(CC[C@@]3(O)C1=CC(=O)[C@]1([H])C[C@@H](O)[C@@H](O)C[C@]21C)[C@@](C)(O)[C@@H](CCC(C)(C)O)OP(O)(O)=O',
                                      'name': '20-hydroxyecdysone 22-phosphate',
                                      'reason': 'No esterified oxygen attached '
                                                'to steroid nucleus found'},
                                  {   'smiles': 'C1=2[C@]3([C@](CC[C@@]1([C@@]4([C@](C(C2)=O)(C[C@H]([C@H](C4)O)O)[H])C)[H])([C@@](CC3)([H])[C@@H]([C@@H](CCC(C)(C)O)OP(=O)([O-])[O-])C)C)O',
                                      'name': 'ecdysone 22-phosphate(2-)',
                                      'reason': 'No esterified oxygen attached '
                                                'to steroid nucleus found'}],
    'attempt': 0,
    'reasoning': 'To classify a molecule as a sterol ester, we need to check '
                 'two main features:\n'
                 '\n'
                 '1. **Presence of a steroid nucleus**: This is characterized '
                 'by a fused ring system consisting of three six-membered '
                 'rings and one five-membered ring (6-6-6-5) in a specific '
                 'fused arrangement.\n'
                 '\n'
                 '2. **Esterification at the 3-hydroxy position**: In sterol '
                 'esters, the hydroxyl group at position 3 of the sterol is '
                 'esterified, forming an ester bond with a fatty acid.\n'
                 '\n'
                 'The following code implements these checks step by step:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 61,
    'num_false_positives': 100,
    'num_true_negatives': 9872,
    'num_false_negatives': 2,
    'num_negatives': None,
    'precision': 0.37888198757763975,
    'recall': 0.9682539682539683,
    'f1': 0.5446428571428571,
    'accuracy': 0.9898355754857997,
    'negative_predictive_value': 0.9997974478428195}