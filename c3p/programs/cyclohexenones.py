"""
Classifies: CHEBI:48953 cyclohexenones
"""
"""
Classifies: cyclohexenones
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cyclohexenones(smiles: str):
    """
    Determines if a molecule is a cyclohexenone based on its SMILES string.
    A cyclohexenone is defined as any six-membered alicyclic ketone having one double bond in the ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclohexenone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all six-membered rings
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    six_membered_rings = [ring for ring in atom_rings if len(ring) == 6]

    if not six_membered_rings:
        return False, "No six-membered rings found"

    # SMARTS pattern for ketone group in ring: [C]=O where C is in ring
    ketone_in_ring_pattern = Chem.MolFromSmarts("[R][C](=O)[R]")

    # Check each six-membered ring
    for ring in six_membered_rings:
        # Check if ring is aliphatic (non-aromatic)
        if any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue  # Skip aromatic rings

        # Get bonds in the ring
        bonds = []
        for i in range(len(ring)):
            bond = mol.GetBondBetweenAtoms(ring[i], ring[(i+1)%len(ring)])
            bonds.append(bond)

        # Count double bonds in the ring
        double_bonds_in_ring = sum(1 for bond in bonds if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE)

        if double_bonds_in_ring != 1:
            continue  # Need exactly one double bond in ring

        # Check for ketone group in the ring
        matches = mol.GetSubstructMatches(ketone_in_ring_pattern)
        for match in matches:
            # Check if both atoms of ketone are in the ring
            if match[0] in ring and match[1] in ring:
                return True, "Contains six-membered alicyclic ring with one double bond and a ketone group in the ring"

    return False, "Does not match criteria for cyclohexenone"

__metadata__ = {
    'chemical_class': {
        'name': 'cyclohexenones',
        'definition': 'Any six-membered alicyclic ketone having one double bond in the ring.'
    }
}


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:48953',
                          'name': 'cyclohexenones',
                          'definition': 'Any six-membered alicyclic ketone '
                                        'having one double bond in the ring.',
                          'parents': ['CHEBI:36132'],
                          'xrefs': ['PMID:21290509'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'CNC(=O)CC[C@H](N)C(O)=O',
                                     'name': 'N(5)-methyl-L-glutamine',
                                     'reason': 'No six-membered rings found'},
                                 {   'smiles': 'C[C@H](OP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c1nc(N)[nH]c2=O)C(O)=O',
                                     'name': "L-lactyl-2-diphospho-5'-guanosine",
                                     'reason': 'Does not match criteria for '
                                               'cyclohexenone'},
                                 {   'smiles': '[H][C@@]1(CC[C@]2(C)[C@]1([H])CC[C@]1([H])[C@@]3(C)CCCC(C)(C)[C@]3([H])CC[C@@]21C)C(=C)CCC=C(C)C',
                                     'name': 'dammara-20,24-diene',
                                     'reason': 'Does not match criteria for '
                                               'cyclohexenone'},
                                 {   'smiles': 'O=C1C2=C(O)C3=C(O)C=CC=C3C(=C2C(C(=O)C)C(C1)(O)C)C=4C=5C(C(O)=C6C4C(C(=O)C)C(O)(C)CC6=O)=C(O)C=CC5',
                                     'name': 'A 39183A',
                                     'reason': 'Does not match criteria for '
                                               'cyclohexenone'},
                                 {   'smiles': 'C[C@H](O)[C@@H](O)C1=Nc2c(NC1)[nH]c(N)nc2=O',
                                     'name': 'L-threo-7,8-dihydrobiopterin',
                                     'reason': 'Does not match criteria for '
                                               'cyclohexenone'},
                                 {   'smiles': 'C[C@@H]1CN([C@H](COC2=C(C=C(C=C2)NC(=O)C)C(=O)N(C[C@H]1OC)C)C)CC3=CC=C(C=C3)F',
                                     'name': 'N-[(4S,7R,8S)-5-[(4-fluorophenyl)methyl]-8-methoxy-4,7,10-trimethyl-11-oxo-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]acetamide',
                                     'reason': 'Does not match criteria for '
                                               'cyclohexenone'},
                                 {   'smiles': 'C(CCN1CCCCC1)(CCN(C(C)C)C(C)=O)(C(N)=O)C2=C(Cl)C=CC=C2',
                                     'name': 'bidisomide',
                                     'reason': 'Does not match criteria for '
                                               'cyclohexenone'},
                                 {   'smiles': 'OCCCCCCCCCCC#CC=C',
                                     'name': '13-Tetradece-11-yn-1-ol',
                                     'reason': 'No six-membered rings found'},
                                 {   'smiles': 'C[C@@H]1CCCCO[C@@H]([C@@H](CN(C(=O)C2=C(O1)C=CC(=C2)N(C)C)[C@@H](C)CO)C)CN(C)CC3=CC=C(C=C3)OC',
                                     'name': '(3R,9S,10R)-16-(dimethylamino)-12-[(2S)-1-hydroxypropan-2-yl]-9-[[(4-methoxyphenyl)methyl-methylamino]methyl]-3,10-dimethyl-2,8-dioxa-12-azabicyclo[12.4.0]octadeca-1(14),15,17-trien-13-one',
                                     'reason': 'Does not match criteria for '
                                               'cyclohexenone'},
                                 {   'smiles': 'O(CCC12CC3CC(C1)CC(C2)C3)C(=O)CC4=CC=C(OCC(O)CNC(C)C)C=C4',
                                     'name': 'Adaprolol',
                                     'reason': 'Does not match criteria for '
                                               'cyclohexenone'}],
    'sample_false_negatives': [],
    'attempt': 0,
    'reasoning': 'To determine if a molecule is a cyclohexenone, we need to '
                 'check for the following:\n'
                 '\n'
                 '- It contains a six-membered alicyclic (non-aromatic) ring.\n'
                 '- The ring has exactly one double bond.\n'
                 '- The ring contains a ketone group (C=O) within the ring.\n'
                 '\n'
                 'We can use RDKit to parse the SMILES string and search for '
                 'these features using SMARTS patterns.\n'
                 '\n'
                 'Here is the code:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 438,
    'num_false_positives': 100,
    'num_true_negatives': 3496,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.8141263940520446,
    'recall': 1.0,
    'f1': 0.8975409836065573,
    'accuracy': 0.9752107089737233,
    'negative_predictive_value': 1.0}