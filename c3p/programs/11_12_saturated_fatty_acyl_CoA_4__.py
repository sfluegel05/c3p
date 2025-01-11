"""
Classifies: CHEBI:84948 11,12-saturated fatty acyl-CoA(4-)
"""
"""
Classifies: 11,12-saturated fatty acyl-CoA(4-)
"""
from rdkit import Chem

def is_11_12_saturated_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 11,12-saturated fatty acyl-CoA(4-) based on its SMILES string.
    Any fatty acyl-CoA(4-) in which the 11-12 bond of the fatty acyl group is saturated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 11,12-saturated fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define thioester pattern to find the fatty acyl-CoA linkage
    thioester_pattern = Chem.MolFromSmarts('[#6](=O)[#16]')  # C(=O)S

    # Find matches of the thioester pattern
    matches = mol.GetSubstructMatches(thioester_pattern)
    if len(matches) == 0:
        return False, "No thioester linkage found (not a fatty acyl-CoA)"

    # Assume the first match is the fatty acyl-CoA linkage
    thioester_carbon_idx = matches[0][0]

    # Function to get the longest carbon chain starting from the carbonyl carbon
    def get_longest_path(mol, start_atom_idx):
        from collections import deque

        visited = set()
        max_path = []

        queue = deque()
        queue.append((start_atom_idx, [start_atom_idx]))

        while queue:
            current_atom_idx, path = queue.popleft()

            atom = mol.GetAtomWithIdx(current_atom_idx)
            neighbors = atom.GetNeighbors()

            terminal_atom = True
            for nbr in neighbors:
                nbr_idx = nbr.GetIdx()
                if nbr_idx in path:
                    continue
                # Proceed only along carbons
                if nbr.GetAtomicNum() != 6:
                    continue
                bond = mol.GetBondBetweenAtoms(current_atom_idx, nbr_idx)
                if bond.GetBondType() not in (Chem.rdchem.BondType.SINGLE, Chem.rdchem.BondType.DOUBLE):
                    continue
                new_path = path + [nbr_idx]
                queue.append((nbr_idx, new_path))
                terminal_atom = False

            if terminal_atom:
                # Reached a terminal carbon
                if len(path) > len(max_path):
                    max_path = path

        return max_path

    # Get the acyl chain as the longest path from the thioester carbon
    acyl_chain = get_longest_path(mol, thioester_carbon_idx)

    if len(acyl_chain) < 12:
        return False, "Fatty acyl chain is too short, less than 12 carbons"

    # Get the atoms at positions 11 and 12 (C11 and C12)
    c11_atom_idx = acyl_chain[10]  # Indexing from 0
    c12_atom_idx = acyl_chain[11]

    # Get the bond between C11 and C12
    bond = mol.GetBondBetweenAtoms(c11_atom_idx, c12_atom_idx)
    if bond is None:
        return False, "No bond between C11 and C12"

    if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
        return True, "C11-C12 bond is saturated (single bond)"
    elif bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
        return False, "C11-C12 bond is unsaturated (double bond)"
    else:
        return False, "C11-C12 bond is neither single nor double"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:84948',
                          'name': '11,12-saturated fatty acyl-CoA(4-)',
                          'definition': 'Any fatty acyl-CoA(4-) in which the '
                                        '11-12 bond of the fatty acyl group is '
                                        'saturated. This is needed for a '
                                        'reaction together with CHEBI:84947',
                          'parents': ['CHEBI:77636'],
                          'xrefs': ['PMID:10563812', 'PMID:2295634'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'CNC(=O)CC[C@H](N)C(O)=O',
                                     'name': 'N(5)-methyl-L-glutamine',
                                     'reason': 'No thioester linkage found '
                                               '(not a fatty acyl-CoA)'},
                                 {   'smiles': 'C[C@H](OP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c1nc(N)[nH]c2=O)C(O)=O',
                                     'name': "L-lactyl-2-diphospho-5'-guanosine",
                                     'reason': 'No thioester linkage found '
                                               '(not a fatty acyl-CoA)'},
                                 {   'smiles': '[H][C@@]1(CC[C@]2(C)[C@]1([H])CC[C@]1([H])[C@@]3(C)CCCC(C)(C)[C@]3([H])CC[C@@]21C)C(=C)CCC=C(C)C',
                                     'name': 'dammara-20,24-diene',
                                     'reason': 'No thioester linkage found '
                                               '(not a fatty acyl-CoA)'},
                                 {   'smiles': 'O=C1C2=C(O)C3=C(O)C=CC=C3C(=C2C(C(=O)C)C(C1)(O)C)C=4C=5C(C(O)=C6C4C(C(=O)C)C(O)(C)CC6=O)=C(O)C=CC5',
                                     'name': 'A 39183A',
                                     'reason': 'No thioester linkage found '
                                               '(not a fatty acyl-CoA)'},
                                 {   'smiles': 'C[C@H](O)[C@@H](O)C1=Nc2c(NC1)[nH]c(N)nc2=O',
                                     'name': 'L-threo-7,8-dihydrobiopterin',
                                     'reason': 'No thioester linkage found '
                                               '(not a fatty acyl-CoA)'},
                                 {   'smiles': 'C[C@@H]1CN([C@H](COC2=C(C=C(C=C2)NC(=O)C)C(=O)N(C[C@H]1OC)C)C)CC3=CC=C(C=C3)F',
                                     'name': 'N-[(4S,7R,8S)-5-[(4-fluorophenyl)methyl]-8-methoxy-4,7,10-trimethyl-11-oxo-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]acetamide',
                                     'reason': 'No thioester linkage found '
                                               '(not a fatty acyl-CoA)'},
                                 {   'smiles': 'C(CCN1CCCCC1)(CCN(C(C)C)C(C)=O)(C(N)=O)C2=C(Cl)C=CC=C2',
                                     'name': 'bidisomide',
                                     'reason': 'No thioester linkage found '
                                               '(not a fatty acyl-CoA)'},
                                 {   'smiles': 'OCCCCCCCCCCC#CC=C',
                                     'name': '13-Tetradece-11-yn-1-ol',
                                     'reason': 'No thioester linkage found '
                                               '(not a fatty acyl-CoA)'},
                                 {   'smiles': 'C[C@@H]1CCCCO[C@@H]([C@@H](CN(C(=O)C2=C(O1)C=CC(=C2)N(C)C)[C@@H](C)CO)C)CN(C)CC3=CC=C(C=C3)OC',
                                     'name': '(3R,9S,10R)-16-(dimethylamino)-12-[(2S)-1-hydroxypropan-2-yl]-9-[[(4-methoxyphenyl)methyl-methylamino]methyl]-3,10-dimethyl-2,8-dioxa-12-azabicyclo[12.4.0]octadeca-1(14),15,17-trien-13-one',
                                     'reason': 'No thioester linkage found '
                                               '(not a fatty acyl-CoA)'},
                                 {   'smiles': 'O(CCC12CC3CC(C1)CC(C2)C3)C(=O)CC4=CC=C(OCC(O)CNC(C)C)C=C4',
                                     'name': 'Adaprolol',
                                     'reason': 'No thioester linkage found '
                                               '(not a fatty acyl-CoA)'}],
    'sample_false_negatives': [],
    'attempt': 0,
    'reasoning': 'To classify molecules as 11,12-saturated fatty acyl-CoA(4-), '
                 'we need to identify the fatty acyl chain attached to the '
                 'Coenzyme A (CoA) moiety via a thioester linkage, and then '
                 'check if the bond between the 11th and 12th carbons in the '
                 'acyl chain is saturated (single bond).\n'
                 '\n'
                 "Here's the step-by-step approach:\n"
                 '\n'
                 '1. **Parse the SMILES string into an RDKit molecule '
                 'object.**\n'
                 '\n'
                 '2. **Identify the thioester linkage** in the molecule using '
                 'a SMARTS pattern that matches a carbonyl carbon attached to '
                 'sulfur (`[C](=O)S`). This carbonyl carbon is the starting '
                 'point (C1) of the fatty acyl chain.\n'
                 '\n'
                 '3. **Traverse the acyl chain** starting from the carbonyl '
                 "carbon, moving along connected carbon atoms. We'll use a "
                 'breadth-first search (BFS) to find the longest carbon chain '
                 'starting from the carbonyl carbon. This helps in cases where '
                 'there might be branching in the molecule.\n'
                 '\n'
                 "4. **Check the length of the acyl chain.** If it's shorter "
                 'than 12 carbons, we cannot check the bond between C11 and '
                 'C12.\n'
                 '\n'
                 '5. **Identify the bond between the 11th and 12th carbons '
                 "(C11 and C12).** We'll get the bond connecting these two "
                 'carbons from our acyl chain path.\n'
                 '\n'
                 "6. **Determine the bond type between C11 and C12.** If it's "
                 'a single bond, the molecule is classified as an '
                 "11,12-saturated fatty acyl-CoA(4-). If it's a double bond, "
                 "it's not.\n"
                 '\n'
                 'Below is the Python code implementing this logic:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 69,
    'num_false_positives': 100,
    'num_true_negatives': 21218,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.40828402366863903,
    'recall': 1.0,
    'f1': 0.5798319327731092,
    'accuracy': 0.995324262402394,
    'negative_predictive_value': 1.0}