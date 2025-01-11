"""
Classifies: CHEBI:65061 2,5-diketopiperazines
"""
"""
Classifies: CHEBI:47909 2,5-diketopiperazine
"""
from rdkit import Chem

def is_2_5_diketopiperazines(smiles: str):
    """
    Determines if a molecule is a 2,5-diketopiperazine based on its SMILES string.
    A 2,5-diketopiperazine has a piperazine-2,5-dione skeleton, which is a six-membered ring 
    with two nitrogen atoms at positions 1 and 4, and two carbonyl groups at positions 2 and 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2,5-diketopiperazine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Iterate over all 6-membered rings
    for ring in atom_rings:
        if len(ring) != 6:
            continue  # Skip rings that are not 6-membered

        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]

        # Count nitrogen atoms in the ring
        num_nitrogens = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 7)

        # Count carbon atoms double-bonded to oxygen (carbonyl groups) in the ring
        num_carbonyl_carbons = 0
        for atom in ring_atoms:
            if atom.GetAtomicNum() == 6:
                # Check if this carbon atom has a double bond to oxygen
                is_carbonyl = False
                for bond in atom.GetBonds():
                    nbr = bond.GetOtherAtom(atom)
                    if nbr.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        is_carbonyl = True
                        break
                if is_carbonyl:
                    num_carbonyl_carbons += 1

        # Check if ring matches 2,5-diketopiperazine skeleton
        if num_nitrogens == 2 and num_carbonyl_carbons == 2:
            return True, "Contains 2,5-diketopiperazine skeleton"

    return False, "Does not contain 2,5-diketopiperazine skeleton"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:47909',
        'name': '2,5-diketopiperazine',
        'definition': 'Any piperazinone that has a piperazine-2,5-dione skeleton.',
        'parents': ['CHEBI:24164', 'CHEBI:48373']
    },
    'config': {
        # Configuration parameters can be added here if needed
    },
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    # Performance metrics can be added here if available
}


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:65061',
                          'name': '2,5-diketopiperazines',
                          'definition': 'Any piperazinone that has a '
                                        'piperazine-2,5-dione skeleton.',
                          'parents': [   'CHEBI:24613',
                                         'CHEBI:46761',
                                         'CHEBI:46846'],
                          'xrefs': ['Wikipedia:Diketopiperazine'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'CNC(=O)CC[C@H](N)C(O)=O',
                                     'name': 'N(5)-methyl-L-glutamine',
                                     'reason': 'Does not contain '
                                               '2,5-diketopiperazine skeleton'},
                                 {   'smiles': 'C[C@H](OP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c1nc(N)[nH]c2=O)C(O)=O',
                                     'name': "L-lactyl-2-diphospho-5'-guanosine",
                                     'reason': 'Does not contain '
                                               '2,5-diketopiperazine skeleton'},
                                 {   'smiles': '[H][C@@]1(CC[C@]2(C)[C@]1([H])CC[C@]1([H])[C@@]3(C)CCCC(C)(C)[C@]3([H])CC[C@@]21C)C(=C)CCC=C(C)C',
                                     'name': 'dammara-20,24-diene',
                                     'reason': 'Does not contain '
                                               '2,5-diketopiperazine skeleton'},
                                 {   'smiles': 'O=C1C2=C(O)C3=C(O)C=CC=C3C(=C2C(C(=O)C)C(C1)(O)C)C=4C=5C(C(O)=C6C4C(C(=O)C)C(O)(C)CC6=O)=C(O)C=CC5',
                                     'name': 'A 39183A',
                                     'reason': 'Does not contain '
                                               '2,5-diketopiperazine skeleton'},
                                 {   'smiles': 'C[C@H](O)[C@@H](O)C1=Nc2c(NC1)[nH]c(N)nc2=O',
                                     'name': 'L-threo-7,8-dihydrobiopterin',
                                     'reason': 'Does not contain '
                                               '2,5-diketopiperazine skeleton'},
                                 {   'smiles': 'C[C@@H]1CN([C@H](COC2=C(C=C(C=C2)NC(=O)C)C(=O)N(C[C@H]1OC)C)C)CC3=CC=C(C=C3)F',
                                     'name': 'N-[(4S,7R,8S)-5-[(4-fluorophenyl)methyl]-8-methoxy-4,7,10-trimethyl-11-oxo-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]acetamide',
                                     'reason': 'Does not contain '
                                               '2,5-diketopiperazine skeleton'},
                                 {   'smiles': 'C(CCN1CCCCC1)(CCN(C(C)C)C(C)=O)(C(N)=O)C2=C(Cl)C=CC=C2',
                                     'name': 'bidisomide',
                                     'reason': 'Does not contain '
                                               '2,5-diketopiperazine skeleton'},
                                 {   'smiles': 'OCCCCCCCCCCC#CC=C',
                                     'name': '13-Tetradece-11-yn-1-ol',
                                     'reason': 'Does not contain '
                                               '2,5-diketopiperazine skeleton'},
                                 {   'smiles': 'C[C@@H]1CCCCO[C@@H]([C@@H](CN(C(=O)C2=C(O1)C=CC(=C2)N(C)C)[C@@H](C)CO)C)CN(C)CC3=CC=C(C=C3)OC',
                                     'name': '(3R,9S,10R)-16-(dimethylamino)-12-[(2S)-1-hydroxypropan-2-yl]-9-[[(4-methoxyphenyl)methyl-methylamino]methyl]-3,10-dimethyl-2,8-dioxa-12-azabicyclo[12.4.0]octadeca-1(14),15,17-trien-13-one',
                                     'reason': 'Does not contain '
                                               '2,5-diketopiperazine skeleton'},
                                 {   'smiles': 'O(CCC12CC3CC(C1)CC(C2)C3)C(=O)CC4=CC=C(OCC(O)CNC(C)C)C=C4',
                                     'name': 'Adaprolol',
                                     'reason': 'Does not contain '
                                               '2,5-diketopiperazine '
                                               'skeleton'}],
    'sample_false_negatives': [],
    'attempt': 0,
    'reasoning': 'To determine if a molecule is a 2,5-diketopiperazine, we '
                 'need to check for the presence of a six-membered ring '
                 'containing two nitrogen atoms and two carbonyl groups (C=O) '
                 'attached to carbons in the ring. We can use RDKit to find '
                 'all six-membered rings in the molecule and test each for '
                 'these features.',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 53,
    'num_false_positives': 100,
    'num_true_negatives': 7345,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.3464052287581699,
    'recall': 1.0,
    'f1': 0.5145631067961165,
    'accuracy': 0.98666311016271,
    'negative_predictive_value': 1.0}