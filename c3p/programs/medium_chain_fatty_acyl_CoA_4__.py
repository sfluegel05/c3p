"""
Classifies: CHEBI:90546 medium-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_medium_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA(4-) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a medium-chain fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the CoA adenylate part and general structure indicative of CoA
    adenylate_pattern = Chem.MolFromSmarts("n1cnc2c1ncnc2N")
    coA_general_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)O")
    if not (mol.HasSubstructMatch(adenylate_pattern) and mol.HasSubstructMatch(coA_general_pattern)):
        return False, "No complete Coenzyme A structure found"

    # Check for the thioester linkage C(=O)S
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if len(thioester_matches) == 0:
        return False, "No thioester bond found"

    # Verify medium-chain fatty acid length (evaluate linear chains starting from carbon after thioester)
    valid_fatty_acid = False
    for match in thioester_matches:
        thioester_carbon_idx = match[1]  # S bound to thioester carbon
        neighbors = mol.GetAtomWithIdx(thioester_carbon_idx).GetNeighbors()
        
        # Find the carbon atom in the acyl chain
        next_c_atom = None
        for neighbor in neighbors:
            if neighbor.GetAtomicNum() == 6:  # Carbon atom
                next_c_atom = neighbor.GetIdx()
                break

        if next_c_atom is not None:
            chain_length = 0
            visited = set()

            # Perform Breadth-First Search (BFS) to traverse carbon chain
            to_visit = [next_c_atom]
            while to_visit:
                atom_idx = to_visit.pop(0)
                if atom_idx in visited:
                    continue
                if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 6:
                    chain_length += 1
                visited.add(atom_idx)
                for neighbor in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
                    if neighbor.GetIdx() not in visited and neighbor.GetAtomicNum() == 6:
                        to_visit.append(neighbor.GetIdx())

            if 6 <= chain_length <= 12:
                valid_fatty_acid = True
                break

    if not valid_fatty_acid:
        return False, "Fatty acid chain not within medium-chain length (6 to 12 carbons)"

    # Confirm the deprotonated state to carry a -4 charge
    total_charge = rdmolops.GetFormalCharge(mol)
    if total_charge != -4:
        return False, f"Total charge is not -4, found {total_charge}"

    return True, "Molecule matches medium-chain fatty acyl-CoA(4-) structure"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:90546',
                          'name': 'medium-chain fatty acyl-CoA(4-)',
                          'definition': 'An acyl-CoA oxoanion that results '
                                        'from deprotonation of the phosphate '
                                        'and diphosphate groups of any '
                                        'medium-chain fatty acyl-CoA; major '
                                        'species at pH 7.3.',
                          'parents': ['CHEBI:77636'],
                          'xrefs': ['MetaCyc:Medium-Chain-Acyl-CoAs'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: NONE\n'
               'False positives: NONE\n'
               'False negatives: SMILES: '
               'CCCCC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12 '
               'NAME: trans-oct-2-enoyl-CoA(4-) REASON: MISSED No complete '
               'Coenzyme A structure found\n'
               ' * SMILES: '
               '[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC([C@H](C(NCCC(NCCSC(=O)/C=C/C=C\\CCCCC)=O)=O)O)(C)C)(=O)[O-])(=O)[O-])[C@H]([C@H]1O)OP([O-])([O-])=O '
               'NAME: (2E,4Z)-deca-2,4-dienoyl-CoA(4-) REASON: MISSED No '
               'complete Coenzyme A structure found\n'
               ' * SMILES: '
               '[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(CC[C@@H](CCCC(C)C)C)=O)=O)=O)O)C)(=O)[O-])(=O)[O-])[C@H]([C@H]1O)OP([O-])([O-])=O '
               'NAME: 4(R),8-dimethylnonanoyl-CoA(4-) REASON: MISSED No '
               'complete Coenzyme A structure found\n'
               ' * SMILES: '
               'CCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12 '
               'NAME: decanoyl-CoA(4-) REASON: MISSED No complete Coenzyme A '
               'structure found\n'
               ' * SMILES: '
               'CCCCCC\\C=C\\CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12 '
               'NAME: trans-dec-3-enoyl-CoA(4-) REASON: MISSED No complete '
               'Coenzyme A structure found\n'
               ' * SMILES: '
               'CCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12 '
               'NAME: heptanoyl-CoA(4-) REASON: MISSED No complete Coenzyme A '
               'structure found\n'
               ' * SMILES: '
               '[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(C(=CCCC)C)=O)=O)=O)O)C)(=O)[O-])(=O)[O-])[C@H]([C@H]1O)OP([O-])([O-])=O '
               'NAME: 2-methylhexenoyl-CoA(4-) REASON: MISSED No complete '
               'Coenzyme A structure found\n'
               ' * SMILES: '
               '[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(CCC(CCCCCC)O)=O)=O)=O)O)C)(=O)[O-])(=O)[O-])[C@H]([C@H]1O)OP([O-])([O-])=O '
               'NAME: 4-hydroxydecanoyl-CoA(4-) REASON: MISSED No complete '
               'Coenzyme A structure found\n'
               ' * SMILES: '
               'CC(C)CCC[C@@H](C)C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12 '
               'NAME: (2R)-2,6-dimethylheptanoyl-CoA(4-) REASON: MISSED No '
               'complete Coenzyme A structure found\n'
               ' * SMILES: '
               '[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(C(CCCCCCCCCC)O)=O)=O)=O)O)C)(=O)[O-])(=O)[O-])[C@H]([C@H]1O)OP([O-])([O-])=O '
               'NAME: 2-hydroxydodecanoyl-CoA(4-) REASON: MISSED No complete '
               'Coenzyme A structure found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'O([C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)CO)[C@H]3[C@@H](O)[C@H](OC(O)[C@@H]3O)CO',
                                     'name': 'beta-D-Glcp-(1->2)-beta-D-Glcp-(1->3)-D-Galp',
                                     'reason': 'No complete Coenzyme A '
                                               'structure found'},
                                 {   'smiles': 'CCS(=O)(=O)NCC[C@@H]1CC[C@@H]([C@H](O1)CO)NC(=O)NC2=CC(=CC(=C2)Cl)Cl',
                                     'name': '1-(3,5-dichlorophenyl)-3-[(2S,3S,6S)-6-[2-(ethylsulfonylamino)ethyl]-2-(hydroxymethyl)-3-oxanyl]urea',
                                     'reason': 'No complete Coenzyme A '
                                               'structure found'},
                                 {   'smiles': 'C(C(O)=O)C/C=C\\C/C=C\\C\\C=C/C=C/C=C/[C@H]1[C@H](C/C=C\\CC)O1',
                                     'name': '(16S,17S)-epoxy-(4Z,7Z,10Z,12E,14E,19Z)-docosahexaenoic '
                                             'acid',
                                     'reason': 'No complete Coenzyme A '
                                               'structure found'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@H]2[C@@H](CO)O[C@@H](O[C@@H]3[C@@H](CO)O[C@@H](O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O',
                                     'name': 'alpha-D-Galp-(1->4)-beta-D-Galp-(1->4)-beta-D-Glcp',
                                     'reason': 'No complete Coenzyme A '
                                               'structure found'},
                                 {   'smiles': '[Li+].[Br-]',
                                     'name': 'lithium bromide',
                                     'reason': 'No complete Coenzyme A '
                                               'structure found'},
                                 {   'smiles': 'C=1(OC)C2=C(C=C3C[C@H]([C@](CC=4C=C(OC)C(OC)=C(C4C13)OC)(C)O)C)OCO2',
                                     'name': 'Besigomsin',
                                     'reason': 'No complete Coenzyme A '
                                               'structure found'},
                                 {   'smiles': 'C1=CC=C(C(=C1)C=CC(=O)C2=CC=CN2)Cl',
                                     'name': '3-(2-chlorophenyl)-1-(1H-pyrrol-2-yl)-2-propen-1-one',
                                     'reason': 'No complete Coenzyme A '
                                               'structure found'},
                                 {   'smiles': 'O=C1OC(O)C(=C1C(O)C(C)C)C',
                                     'name': '5-Hydroxy-3-(1-hydroxy-2-methylpropyl)-4-methyl-2(5H)-furanone',
                                     'reason': 'No complete Coenzyme A '
                                               'structure found'},
                                 {   'smiles': '[H][C@@]1(COC(C)(C)[C@@]1([H])OC(=O)\\C=C\\c1ccccc1)c1c(O)ccc2C(=O)C[C@H](Oc12)c1ccccc1',
                                     'name': '(+)-tephrorin B',
                                     'reason': 'No complete Coenzyme A '
                                               'structure found'},
                                 {   'smiles': 'C[C@@H]1CN(C(=O)C2=C(C3=CC=CC=C3CO[C@@H]1CN(C)CC4=CC=CC(=C4)C(=O)O)C5=CC=CC=C5N2C)[C@@H](C)CO',
                                     'name': 'LSM-9341',
                                     'reason': 'No complete Coenzyme A '
                                               'structure found'}],
    'sample_false_negatives': [   {   'smiles': 'C(C(CCCCCCCCC)C)C(SCCNC(CCNC(=O)[C@@H](C(COP(OP(OC[C@H]1O[C@@H](N2C3=C(C(=NC=N3)N)N=C2)[C@@H]([C@@H]1OP([O-])([O-])=O)O)(=O)[O-])(=O)[O-])(C)C)O)=O)=O',
                                      'name': '3-methyldodecanoyl-CoA(4-)',
                                      'reason': 'Fatty acid chain not within '
                                                'medium-chain length (6 to 12 '
                                                'carbons)'},
                                  {   'smiles': '[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(CCC(OP([O-])([O-])=O)CCCCCC)=O)=O)=O)O)C)(=O)[O-])(=O)[O-])[C@H]([C@H]1O)OP([O-])([O-])=O',
                                      'name': '4-phosphodecanoyl-CoA(6-)',
                                      'reason': 'Total charge is not -4, found '
                                                '-6'}],
    'attempt': 3,
    'reasoning': 'To create a program that correctly classifies medium-chain '
                 'fatty acyl-CoA(4-) molecules based on their SMILES strings, '
                 'we need to carefully examine the definition and '
                 'characteristics of this chemical class. A medium-chain fatty '
                 'acyl-CoA(4-) is defined as an acyl-CoA oxoanion with '
                 'deprotonated phosphate and diphosphate groups, typically '
                 'characterized by fatty acyl chains with 6-12 carbon atoms.\n'
                 '\n'
                 "Let's analyze the previous program and the outcomes to "
                 'identify potential improvements:\n'
                 '\n'
                 '### Issues & Hypotheses:\n'
                 '\n'
                 '1. **CoA Structure Identification:**\n'
                 '   - The previous program may not accurately capture the '
                 'complex Coenzyme A structure. We must ensure the SMARTS '
                 'pattern represents the full CoA structure, which includes '
                 "adenosine 3'-phosphate 5'-diphosphate linked to pantetheine "
                 'with a thioester linkage to a fatty acyl group.\n'
                 '\n'
                 '2. **Medium-Chain Definition:**\n'
                 '   - Ensure the program verifies the fatty acyl chain length '
                 '(6-12 carbons) from the thioester carbonyl correctly. The '
                 'method of iterating from the thioester may not be robust '
                 'enough and could fail to accurately count the carbon atoms, '
                 'especially if structures have branching or unsaturation '
                 'within this range.\n'
                 '\n'
                 '3. **Charge Verification:**\n'
                 '   - The program should properly calculate and check the '
                 "molecule's formal charge, ensuring deprotonation accounts "
                 'for the -4 charge.\n'
                 '\n'
                 '### Improvement Plan:\n'
                 '\n'
                 '1. Use a more robust way to recognize the entire CoA '
                 'structure, particularly focusing on the adenine, ribose, '
                 'diphosphate, pantetheine, and thioester connections.\n'
                 '2. Guarantee the chain length calculation from the thioester '
                 'carbonyl identifies the complete linear or branched chain '
                 'accurately.\n'
                 '3. Validate the formal charge precisely and ensure '
                 'deprotonation of the appropriate groups.\n'
                 '4. Add detailed logging for debugging when matches fail, '
                 'especially for complex biochemical structures like CoA.\n'
                 '\n'
                 '### Updated Code:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 43,
    'num_false_positives': 100,
    'num_true_negatives': 102281,
    'num_false_negatives': 2,
    'num_negatives': None,
    'precision': 0.3006993006993007,
    'recall': 0.9555555555555556,
    'f1': 0.45744680851063824,
    'accuracy': 0.9990041591002284,
    'negative_predictive_value': 0.999980446408494}