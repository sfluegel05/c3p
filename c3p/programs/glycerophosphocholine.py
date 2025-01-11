"""
Classifies: CHEBI:36313 glycerophosphocholine
"""
"""
Classifies: glycerophosphocholine
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glycerophosphocholine(smiles: str):
    """
    Determines if a molecule is a glycerophosphocholine based on its SMILES string.
    A glycerophosphocholine has a glycerol backbone with fatty acid chains attached
    at positions sn-1 and sn-2 via ester or ether bonds, and a phosphocholine group at position sn-3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophosphocholine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns
    # Glycerol backbone with positions sn-1, sn-2, sn-3
    glycerol_pattern = Chem.MolFromSmarts("[O][CH2][CH](O)[CH2][O]")
    if glycerol_pattern is None:
        return False, "Invalid glycerol backbone SMARTS pattern"

    # Phosphocholine group
    phosphocholine_pattern = Chem.MolFromSmarts("O[P](=O)([O-])OCC[N+](C)(C)C")
    if phosphocholine_pattern is None:
        return False, "Invalid phosphocholine SMARTS pattern"

    # Check for glycerol backbone
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not glycerol_matches:
        return False, "No glycerol backbone found"

    # Check for phosphocholine group
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine group found"

    # Check for fatty acid chains attached via ester or ether bonds at sn-1 and sn-2
    # Ester or ether bonds can be represented by O-C pattern
    ester_ether_pattern = Chem.MolFromSmarts("O[C]")
    if ester_ether_pattern is None:
        return False, "Invalid ester/ether SMARTS pattern"

    fatty_acid_chains = 0
    # Check attachments at sn-1 and sn-2 positions in the glycerol backbone
    for match in glycerol_matches:
        sn1_O_idx = match[0]  # Oxygen at sn-1
        sn2_O_idx = match[2]  # Oxygen at sn-2

        # Check for ester or ether at sn-1
        sn1_O_atom = mol.GetAtomWithIdx(sn1_O_idx)
        sn1_attached = False
        for bond in sn1_O_atom.GetBonds():
            neighbor_atom = bond.GetOtherAtom(sn1_O_atom)
            if neighbor_atom.GetIdx() != match[1]:  # Exclude bond to glycerol backbone
                if mol.GetSubstructMatch(ester_ether_pattern, useChirality=False):
                    fatty_acid_chains += 1
                    sn1_attached = True
                    break

        # Check for ester or ether at sn-2
        sn2_O_atom = mol.GetAtomWithIdx(sn2_O_idx)
        sn2_attached = False
        for bond in sn2_O_atom.GetBonds():
            neighbor_atom = bond.GetOtherAtom(sn2_O_atom)
            if neighbor_atom.GetIdx() != match[1]:  # Exclude bond to glycerol backbone
                if mol.GetSubstructMatch(ester_ether_pattern, useChirality=False):
                    fatty_acid_chains += 1
                    sn2_attached = True
                    break

        if fatty_acid_chains >= 2:
            return True, "Molecule is a glycerophosphocholine"

    return False, "Does not match glycerophosphocholine structure"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'glycerophosphocholine',
        'definition': 'The glycerol phosphate ester of a phosphocholine. A nutrient with many different roles in human health.',
        'parents': []
    }
}


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36313',
                          'name': 'glycerophosphocholine',
                          'definition': 'The glycerol phosphate ester of a '
                                        'phosphocholine. A nutrient with many '
                                        'different roles in human health.',
                          'parents': ['CHEBI:36700', 'CHEBI:37739'],
                          'xrefs': ['PMID:8467564'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Error: Python argument types in\n'
               '    Mol.GetSubstructMatches(Mol, NoneType)\n'
               'did not match C++ signature:\n'
               '    GetSubstructMatches(RDKit::ROMol self, RDKit::MolBundle '
               'query, RDKit::SubstructMatchParameters params)\n'
               '    GetSubstructMatches(RDKit::ROMol self, RDKit::ROMol query, '
               'RDKit::SubstructMatchParameters params)\n'
               '    GetSubstructMatches(RDKit::ROMol self, RDKit::MolBundle '
               'query, bool uniquify=True, bool useChirality=False, bool '
               'useQueryQueryMatches=False, unsigned int maxMatches=1000)\n'
               '    GetSubstructMatches(RDKit::ROMol self, RDKit::ROMol query, '
               'bool uniquify=True, bool useChirality=False, bool '
               'useQueryQueryMatches=False, unsigned int maxMatches=1000)\n'
               'Attempt failed: F1 score of 0 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: NONE\n'
               'False positives: NONE\n'
               'False negatives: NONE\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'CNC(=O)CC[C@H](N)C(O)=O',
                                     'name': 'N(5)-methyl-L-glutamine',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'C[C@H](OP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c1nc(N)[nH]c2=O)C(O)=O',
                                     'name': "L-lactyl-2-diphospho-5'-guanosine",
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': '[H][C@@]1(CC[C@]2(C)[C@]1([H])CC[C@]1([H])[C@@]3(C)CCCC(C)(C)[C@]3([H])CC[C@@]21C)C(=C)CCC=C(C)C',
                                     'name': 'dammara-20,24-diene',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'O=C1C2=C(O)C3=C(O)C=CC=C3C(=C2C(C(=O)C)C(C1)(O)C)C=4C=5C(C(O)=C6C4C(C(=O)C)C(O)(C)CC6=O)=C(O)C=CC5',
                                     'name': 'A 39183A',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'C[C@H](O)[C@@H](O)C1=Nc2c(NC1)[nH]c(N)nc2=O',
                                     'name': 'L-threo-7,8-dihydrobiopterin',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'C[C@@H]1CN([C@H](COC2=C(C=C(C=C2)NC(=O)C)C(=O)N(C[C@H]1OC)C)C)CC3=CC=C(C=C3)F',
                                     'name': 'N-[(4S,7R,8S)-5-[(4-fluorophenyl)methyl]-8-methoxy-4,7,10-trimethyl-11-oxo-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]acetamide',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'C(CCN1CCCCC1)(CCN(C(C)C)C(C)=O)(C(N)=O)C2=C(Cl)C=CC=C2',
                                     'name': 'bidisomide',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'OCCCCCCCCCCC#CC=C',
                                     'name': '13-Tetradece-11-yn-1-ol',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'C[C@@H]1CCCCO[C@@H]([C@@H](CN(C(=O)C2=C(O1)C=CC(=C2)N(C)C)[C@@H](C)CO)C)CN(C)CC3=CC=C(C=C3)OC',
                                     'name': '(3R,9S,10R)-16-(dimethylamino)-12-[(2S)-1-hydroxypropan-2-yl]-9-[[(4-methoxyphenyl)methyl-methylamino]methyl]-3,10-dimethyl-2,8-dioxa-12-azabicyclo[12.4.0]octadeca-1(14),15,17-trien-13-one',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'O(CCC12CC3CC(C1)CC(C2)C3)C(=O)CC4=CC=C(OCC(O)CNC(C)C)C=C4',
                                     'name': 'Adaprolol',
                                     'reason': 'No glycerol backbone found'}],
    'sample_false_negatives': [   {   'smiles': 'CCCCCCCCCCCCCCCCSC[C@H](COP(O)(=O)OCC[N+](C)(C)C)OP(O)(=O)CCCCCCCCCCCCCCCC',
                                      'name': '1-S-hexadecyl-2-O-[hexadecyl(hydroxy)phosphoryl]-1-thio-sn-glycero-3-phosphocholine',
                                      'reason': 'No glycerol backbone found'},
                                  {   'smiles': 'CCCCCCCCCCCC\\C=C/OC[C@H](COP(O)(=O)OCC[N+](C)(C)C)OC(=O)CCCCCCC\\C=C/CCCCCCCC',
                                      'name': '1-O-[(Z)-tetradec-1-enyl]-2-O-[(Z)-octadec-9-enoyl]-sn-glycero-3-phosphocholine',
                                      'reason': 'No phosphocholine group '
                                                'found'},
                                  {   'smiles': 'P(OCC[N+](C)(C)C)(OC[C@H](O)COC(=O)CC=CCC=CCCCCCCCCCCCCCCCCC)(O)=O',
                                      'name': 'LPC(24:2)',
                                      'reason': 'No phosphocholine group '
                                                'found'},
                                  {   'smiles': 'OC[C@H](COP(OCC[N+](C)(C)C)(=O)[O-])OC(=O)CCCCCCC/C=C\\CCCCCC',
                                      'name': '2-[(9Z)-hexadecenoyl]-sn-glycero-3-phosphocholine',
                                      'reason': 'Does not match '
                                                'glycerophosphocholine '
                                                'structure'},
                                  {   'smiles': 'C(CCCCCCCCCCCC)CCCSC[C@@H](NC(CCCCCCCCCCCCCCC)=O)COP(OCC[N+](C)(C)C)(=O)[O-]',
                                      'name': 'thioetheramide PC',
                                      'reason': 'No glycerol backbone found'},
                                  {   'smiles': 'P(OCC[N+](C)(C)C)(OC[C@H](O)COCCCCCCCCCCCCCCCC)(O)=O',
                                      'name': 'Lyso-PAF C-16-d4',
                                      'reason': 'No phosphocholine group '
                                                'found'},
                                  {   'smiles': 'P(OCC[N+](C)(C)C)(OCC(OC(=O)CCCCCCCCC=1OC(=C(C1C)C)CCC)COC(=O)CCCCCCCCC=2OC(=C(C2C)C)CCCCC)(O)=O',
                                      'name': '(3-{[9-(3,4-dimethyl-5-pentylfuran-2-yl)nonanoyl]oxy}-2-{[9-(3,4-dimethyl-5-propylfuran-2-yl)nonanoyl]oxy}propoxy)[2-(trimethylazaniumyl)ethoxy]phosphinic '
                                              'acid',
                                      'reason': 'No phosphocholine group '
                                                'found'},
                                  {   'smiles': 'P(OCC[N+](C)(C)C)(OC[C@H](O)CO)(O)=O',
                                      'name': 'Glycero-3-phosphocholine',
                                      'reason': 'No phosphocholine group '
                                                'found'},
                                  {   'smiles': 'P(OC[C@H](OC(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCCCCCCCC/C=C\\C/C=C\\CCCCC)(O)(O)=O',
                                      'name': 'PA(20:2(11Z,14Z)/22:6(4Z,7Z,10Z,13Z,16Z,19Z))',
                                      'reason': 'No phosphocholine group '
                                                'found'},
                                  {   'smiles': 'P(OCC[N+](C)(C)C)(OC[C@H](O)COC(=O)CC=CCC=CCC=CCC=CCCCCCCCCCCC)(O)=O',
                                      'name': 'LPC(24:4)',
                                      'reason': 'No phosphocholine group '
                                                'found'}],
    'attempt': 2,
    'reasoning': 'The error in the previous program was caused by invalid '
                 'SMARTS patterns leading to `Chem.MolFromSmarts` returning '
                 '`None`. Specifically, using triple-quoted multi-line strings '
                 'for SMARTS patterns can introduce unintended whitespace or '
                 'line breaks, which makes the SMARTS invalid. Additionally, '
                 'the SMARTS syntax may have been incorrect.\n'
                 '\n'
                 "To fix this, we'll:\n"
                 '\n'
                 '- Use single-line strings for SMARTS patterns without '
                 'unnecessary whitespace or newlines.\n'
                 '- Ensure the SMARTS syntax is correct.\n'
                 '- Add checks to verify that `Chem.MolFromSmarts` does not '
                 'return `None`.\n'
                 '- Simplify the approach by defining clear patterns for the '
                 'glycerol backbone, fatty acid chains, and phosphocholine '
                 'group.\n'
                 '\n'
                 'Here is the corrected code:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 979,
    'num_false_positives': 30,
    'num_true_negatives': 141256,
    'num_false_negatives': 35,
    'num_negatives': None,
    'precision': 0.9702675916749257,
    'recall': 0.965483234714004,
    'f1': 0.967869500741473,
    'accuracy': 0.9995432185523542,
    'negative_predictive_value': 0.9997522842927009}