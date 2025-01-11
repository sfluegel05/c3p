"""
Classifies: CHEBI:61907 medium-chain fatty acyl-CoA
"""
"""
Classifies: medium-chain fatty acyl-CoA
Definition: A fatty acyl-CoA that results from the formal condensation of the thiol group of coenzyme A with the carboxy group of any medium-chain fatty acid.
"""

from rdkit import Chem

def is_medium_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a medium-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify the thioester linkage: C(=O)S
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if len(thioester_matches) != 1:
        return False, f"Expected one thioester linkage, found {len(thioester_matches)}"

    # Get indices of carbonyl carbon and sulfur atom in the thioester linkage
    carbonyl_c_idx = thioester_matches[0][0]
    sulfur_idx = thioester_matches[0][2]

    # Collect atoms of the acyl chain excluding the sulfur atom
    visited_atoms = set()

    def traverse_acyl_chain(atom_idx, exclude_idx):
        """
        Recursively traverse the acyl chain starting from the carbonyl carbon,
        excluding the path towards the sulfur atom (CoA moiety).
        """
        if atom_idx in visited_atoms:
            return
        visited_atoms.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx != exclude_idx:
                traverse_acyl_chain(neighbor_idx, exclude_idx)

    # Start traversal from the carbonyl carbon, exclude sulfur atom
    traverse_acyl_chain(carbonyl_c_idx, sulfur_idx)

    # Count the number of carbons in the acyl chain
    num_carbons = sum(1 for idx in visited_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)

    # Check if the number of carbons is within medium-chain length (6-12)
    if num_carbons < 6:
        return False, f"Acyl chain too short ({num_carbons} carbons), not medium-chain"
    if num_carbons > 12:
        return False, f"Acyl chain too long ({num_carbons} carbons), not medium-chain"

    # Verify the presence of coenzyme A moiety
    # Define a SMARTS pattern for key features of CoA
    coa_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(O)=O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A moiety not found"

    return True, f"Contains medium-chain acyl group ({num_carbons} carbons) attached to coenzyme A via thioester linkage"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:61907',
                          'name': 'medium-chain fatty acyl-CoA',
                          'definition': 'A fatty acyl-CoA that results from '
                                        'the formal condensation of the thiol '
                                        'group of coenzyme A with the carboxy '
                                        'group of any medium-chain fatty acid.',
                          'parents': ['CHEBI:37554'],
                          'xrefs': [   'PMID:11027146',
                                       'PMID:1309181',
                                       'PMID:19170545',
                                       'PMID:8218225',
                                       'PMID:8267794'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'CNC(=O)CC[C@H](N)C(O)=O',
                                     'name': 'N(5)-methyl-L-glutamine',
                                     'reason': 'Expected one thioester '
                                               'linkage, found 0'},
                                 {   'smiles': 'C[C@H](OP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c1nc(N)[nH]c2=O)C(O)=O',
                                     'name': "L-lactyl-2-diphospho-5'-guanosine",
                                     'reason': 'Expected one thioester '
                                               'linkage, found 0'},
                                 {   'smiles': '[H][C@@]1(CC[C@]2(C)[C@]1([H])CC[C@]1([H])[C@@]3(C)CCCC(C)(C)[C@]3([H])CC[C@@]21C)C(=C)CCC=C(C)C',
                                     'name': 'dammara-20,24-diene',
                                     'reason': 'Expected one thioester '
                                               'linkage, found 0'},
                                 {   'smiles': 'O=C1C2=C(O)C3=C(O)C=CC=C3C(=C2C(C(=O)C)C(C1)(O)C)C=4C=5C(C(O)=C6C4C(C(=O)C)C(O)(C)CC6=O)=C(O)C=CC5',
                                     'name': 'A 39183A',
                                     'reason': 'Expected one thioester '
                                               'linkage, found 0'},
                                 {   'smiles': 'C[C@H](O)[C@@H](O)C1=Nc2c(NC1)[nH]c(N)nc2=O',
                                     'name': 'L-threo-7,8-dihydrobiopterin',
                                     'reason': 'Expected one thioester '
                                               'linkage, found 0'},
                                 {   'smiles': 'C[C@@H]1CN([C@H](COC2=C(C=C(C=C2)NC(=O)C)C(=O)N(C[C@H]1OC)C)C)CC3=CC=C(C=C3)F',
                                     'name': 'N-[(4S,7R,8S)-5-[(4-fluorophenyl)methyl]-8-methoxy-4,7,10-trimethyl-11-oxo-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]acetamide',
                                     'reason': 'Expected one thioester '
                                               'linkage, found 0'},
                                 {   'smiles': 'C(CCN1CCCCC1)(CCN(C(C)C)C(C)=O)(C(N)=O)C2=C(Cl)C=CC=C2',
                                     'name': 'bidisomide',
                                     'reason': 'Expected one thioester '
                                               'linkage, found 0'},
                                 {   'smiles': 'OCCCCCCCCCCC#CC=C',
                                     'name': '13-Tetradece-11-yn-1-ol',
                                     'reason': 'Expected one thioester '
                                               'linkage, found 0'},
                                 {   'smiles': 'C[C@@H]1CCCCO[C@@H]([C@@H](CN(C(=O)C2=C(O1)C=CC(=C2)N(C)C)[C@@H](C)CO)C)CN(C)CC3=CC=C(C=C3)OC',
                                     'name': '(3R,9S,10R)-16-(dimethylamino)-12-[(2S)-1-hydroxypropan-2-yl]-9-[[(4-methoxyphenyl)methyl-methylamino]methyl]-3,10-dimethyl-2,8-dioxa-12-azabicyclo[12.4.0]octadeca-1(14),15,17-trien-13-one',
                                     'reason': 'Expected one thioester '
                                               'linkage, found 0'},
                                 {   'smiles': 'O(CCC12CC3CC(C1)CC(C2)C3)C(=O)CC4=CC=C(OCC(O)CNC(C)C)C=C4',
                                     'name': 'Adaprolol',
                                     'reason': 'Expected one thioester '
                                               'linkage, found 0'}],
    'sample_false_negatives': [   {   'smiles': 'C(C(CCCCCCCCC)C)C(SCCNC(CCNC(=O)[C@@H](C(COP(OP(OC[C@H]1O[C@@H](N2C3=C(C(=NC=N3)N)N=C2)[C@@H]([C@@H]1OP(O)(O)=O)O)(=O)O)(=O)O)(C)C)O)=O)=O',
                                      'name': '3-methyldodecanoyl-CoA',
                                      'reason': 'Acyl chain too long (13 '
                                                'carbons), not medium-chain'},
                                  {   'smiles': 'CC(C)CCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12',
                                      'name': 'isotridecanoyl-CoA',
                                      'reason': 'Acyl chain too long (13 '
                                                'carbons), not medium-chain'}],
    'attempt': 0,
    'reasoning': 'To determine if a molecule is a medium-chain fatty acyl-CoA, '
                 'we need to:\n'
                 '\n'
                 '1. **Verify the presence of a fatty acyl-CoA moiety**: The '
                 'molecule should contain a coenzyme A (CoA) attached via a '
                 'thioester linkage to an acyl group.\n'
                 '\n'
                 '2. **Determine the acyl chain length**: The acyl chain '
                 'attached to CoA should be derived from a medium-chain fatty '
                 'acid, typically with 6 to 12 carbons.\n'
                 '\n'
                 'We can implement this by:\n'
                 '\n'
                 '- Parsing the SMILES string to create a molecule object.\n'
                 '- Identifying the thioester linkage (`C(=O)S`) in the '
                 'molecule.\n'
                 '- Traversing the acyl chain starting from the carbonyl '
                 'carbon, excluding the sulfur atom (which leads to CoA), to '
                 'collect all connected carbon atoms.\n'
                 '- Counting the number of carbons in the acyl chain.\n'
                 '- Checking if the acyl chain length falls within the '
                 'medium-chain range (6-12 carbons).\n'
                 '- Verifying the presence of the CoA moiety by matching key '
                 'structural features.',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 46,
    'num_false_positives': 100,
    'num_true_negatives': 33838,
    'num_false_negatives': 2,
    'num_negatives': None,
    'precision': 0.3150684931506849,
    'recall': 0.9583333333333334,
    'f1': 0.4742268041237113,
    'accuracy': 0.9969987641970223,
    'negative_predictive_value': 0.9999408983451537}