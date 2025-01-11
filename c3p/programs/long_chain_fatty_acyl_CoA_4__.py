"""
Classifies: CHEBI:83139 long-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem

def is_long_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA(4-) based on its SMILES string.
    The compound should have a long-chain fatty acid linked to Coenzyme A via a thioester bond,
    with deprotonated phosphate groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the CoA thioester linkage pattern
    # The pattern for the thioester linkage of CoA to the fatty acid chain
    thioester_patterns = [
        Chem.MolFromSmarts("C(=O)SCCNC(=O)CCNC(=O)C"),
        Chem.MolFromSmarts("C(=O)SCCNC(=O)CCN")
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in thioester_patterns):
        return False, "No CoA thioester linkage found"

    # Check for a long fatty acid chain, allowing for both straight and unsaturated chains
    long_chain_patterns = [
        Chem.MolFromSmarts("C(=O)SCCC[C;R0]"),  # Long chain allowing flexibility
        Chem.MolFromSmarts("[C;R0]C(=O)SCCNC(=O)CCN")  # Include variations in long chain lengths
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in long_chain_patterns):
        return False, "Insufficient carbon chain length or structure in the fatty acid"

    # Look for the CoA backbone with deprotonated phosphate groups
    coa_diphosphate_pattern = Chem.MolFromSmarts("COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)")
    if not mol.HasSubstructMatch(coa_diphosphate_pattern):
        return False, "CoA backbone or deprotonated phosphate groups not found"

    return True, "Matches long-chain fatty acyl-CoA(4-) structure with deprotonated phosphate groups"

# Example usage
smiles_example = "CC(C)(COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12)[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)CCCCCCC\C=C/C\C=C/CC=C"
result, reason = is_long_chain_fatty_acyl_CoA_4__(smiles_example)
print(result, reason)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:83139',
                          'name': 'long-chain fatty acyl-CoA(4-)',
                          'definition': 'A fatty acyl-CoA(4-) arising from '
                                        'deprotonation of the phosphate and '
                                        'diphosphate OH groups of any '
                                        'long-chain fatty acyl-CoA; major '
                                        'species at pH 7.3.',
                          'parents': ['CHEBI:77636'],
                          'xrefs': ['MetaCyc:Long-Chain-Acyl-CoAs'],
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
               'CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/C=C/CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12 '
               'NAME: (3E,5Z,8Z,11Z,14Z)-icosapentaenoyl-CoA(4-) REASON: '
               'MISSED Insufficient carbon chain length in the fatty acid\n'
               ' * SMILES: '
               'CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12 '
               'NAME: (4Z,7Z,10Z,13Z,16Z,19Z)-docosahexaenoyl-CoA(4-) REASON: '
               'MISSED Insufficient carbon chain length in the fatty acid\n'
               ' * SMILES: '
               'CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12 '
               'NAME: (3R,7Z,10Z,13Z,16Z)-3-hydroxydocosatetraenoyl-CoA(4-) '
               'REASON: MISSED Insufficient carbon chain length in the fatty '
               'acid\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12 '
               'NAME: pentadecanoyl-CoA(4-) REASON: MISSED Insufficient carbon '
               'chain length in the fatty acid\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCC(O)C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12 '
               'NAME: 2-hydroxypalmitoyl-CoA(4-) REASON: MISSED Insufficient '
               'carbon chain length in the fatty acid\n'
               ' * SMILES: '
               'CCCCCCCCCCCC(O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12 '
               'NAME: 3-hydroxytetradecanoyl-CoA(4-) REASON: MISSED '
               'Insufficient carbon chain length in the fatty acid\n'
               ' * SMILES: '
               'CCCCC\\C=C/C\\C=C/C\\C=C/CCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12 '
               'NAME: (10Z,13Z,16Z)-docosatrienoyl-CoA(4-) REASON: MISSED '
               'Insufficient carbon chain length in the fatty acid\n'
               ' * SMILES: '
               'CCCCCC\\C=C\\C=C/CCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12 '
               'NAME: (9Z,11E)-octadecadienoyl-CoA(4-) REASON: MISSED '
               'Insufficient carbon chain length in the fatty acid\n'
               ' * SMILES: '
               'CCCCC\\C=C/CCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12 '
               'NAME: (9Z)-pentadecenoyl-CoA(4-) REASON: MISSED Insufficient '
               'carbon chain length in the fatty acid\n'
               ' * SMILES: '
               'CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12 '
               'NAME: (7Z,10Z,13Z,16Z)-3-oxodocosatetraenoyl-CoA(4-) REASON: '
               'MISSED Insufficient carbon chain length in the fatty acid\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'O([C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)CO)[C@H]3[C@@H](O)[C@H](OC(O)[C@@H]3O)CO',
                                     'name': 'beta-D-Glcp-(1->2)-beta-D-Glcp-(1->3)-D-Galp',
                                     'reason': 'No CoA thioester linkage '
                                               'found'},
                                 {   'smiles': 'CCS(=O)(=O)NCC[C@@H]1CC[C@@H]([C@H](O1)CO)NC(=O)NC2=CC(=CC(=C2)Cl)Cl',
                                     'name': '1-(3,5-dichlorophenyl)-3-[(2S,3S,6S)-6-[2-(ethylsulfonylamino)ethyl]-2-(hydroxymethyl)-3-oxanyl]urea',
                                     'reason': 'No CoA thioester linkage '
                                               'found'},
                                 {   'smiles': 'C(C(O)=O)C/C=C\\C/C=C\\C\\C=C/C=C/C=C/[C@H]1[C@H](C/C=C\\CC)O1',
                                     'name': '(16S,17S)-epoxy-(4Z,7Z,10Z,12E,14E,19Z)-docosahexaenoic '
                                             'acid',
                                     'reason': 'No CoA thioester linkage '
                                               'found'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@H]2[C@@H](CO)O[C@@H](O[C@@H]3[C@@H](CO)O[C@@H](O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O',
                                     'name': 'alpha-D-Galp-(1->4)-beta-D-Galp-(1->4)-beta-D-Glcp',
                                     'reason': 'No CoA thioester linkage '
                                               'found'},
                                 {   'smiles': '[Li+].[Br-]',
                                     'name': 'lithium bromide',
                                     'reason': 'No CoA thioester linkage '
                                               'found'},
                                 {   'smiles': 'C=1(OC)C2=C(C=C3C[C@H]([C@](CC=4C=C(OC)C(OC)=C(C4C13)OC)(C)O)C)OCO2',
                                     'name': 'Besigomsin',
                                     'reason': 'No CoA thioester linkage '
                                               'found'},
                                 {   'smiles': 'C1=CC=C(C(=C1)C=CC(=O)C2=CC=CN2)Cl',
                                     'name': '3-(2-chlorophenyl)-1-(1H-pyrrol-2-yl)-2-propen-1-one',
                                     'reason': 'No CoA thioester linkage '
                                               'found'},
                                 {   'smiles': 'O=C1OC(O)C(=C1C(O)C(C)C)C',
                                     'name': '5-Hydroxy-3-(1-hydroxy-2-methylpropyl)-4-methyl-2(5H)-furanone',
                                     'reason': 'No CoA thioester linkage '
                                               'found'},
                                 {   'smiles': '[H][C@@]1(COC(C)(C)[C@@]1([H])OC(=O)\\C=C\\c1ccccc1)c1c(O)ccc2C(=O)C[C@H](Oc12)c1ccccc1',
                                     'name': '(+)-tephrorin B',
                                     'reason': 'No CoA thioester linkage '
                                               'found'},
                                 {   'smiles': 'C[C@@H]1CN(C(=O)C2=C(C3=CC=CC=C3CO[C@@H]1CN(C)CC4=CC=CC(=C4)C(=O)O)C5=CC=CC=C5N2C)[C@@H](C)CO',
                                     'name': 'LSM-9341',
                                     'reason': 'No CoA thioester linkage '
                                               'found'}],
    'sample_false_negatives': [   {   'smiles': '[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC([C@H](C(NCCC(NCCSC(=O)C[C@@H](O)CCC/C=C\\CCCCCC)=O)=O)O)(C)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O',
                                      'name': '(3S,7Z)-3-hydroxytetradec-7-enoyl-CoA(4-)',
                                      'reason': 'CoA backbone or deprotonated '
                                                'phosphate groups not found'}],
    'attempt': 4,
    'reasoning': 'To address the issues identified in the previous '
                 "implementation, let's analyze and refine our approach for "
                 'detecting long-chain fatty acyl-CoA(4-) compounds based on '
                 "their SMILES representation. Here's a breakdown of what was "
                 'likely missed and how we can improve:\n'
                 '\n'
                 '### Problem Analysis\n'
                 '1. **Chain Length Detection**: The previous implementation '
                 'checked for a long chain of at least 12 carbons but failed '
                 'to recognize valid compounds due to incorrect SMILES pattern '
                 'recognition or checking the wrong entry point of the chain.\n'
                 '\n'
                 '2. **Phosphate Group Detection**: The SMARTS pattern for '
                 'recognizing the deprotonated phosphate groups was possibly '
                 'too stringent or incorrect, leading to missed detections. We '
                 "need to ensure we're capturing a broad and adaptable "
                 "representation that aligns with the compounds' common "
                 'patterns.\n'
                 '\n'
                 '3. **Flexibility and Variability**: Long-chain fatty '
                 'acyl-CoAs might have variations in their structures, such as '
                 'unsaturated bonds or specific branching, not covered in the '
                 'initial implementation.\n'
                 '\n'
                 '### Solutions\n'
                 '1. **Flexible SMARTS**: Allow some variability in the '
                 'patterns for carbon chains and identify ester or thioester '
                 'patterns more effectively.\n'
                 '   \n'
                 '2. **Review Patterns**: Ensure the SMARTS patterns for both '
                 'the CoA linkage and phosphate groups capture the typical '
                 'core structure while allowing for variations in the chain '
                 'length and unsaturation levels.\n'
                 '\n'
                 "Here's the revised implementation with these considerations:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 136,
    'num_false_positives': 100,
    'num_true_negatives': 29749,
    'num_false_negatives': 1,
    'num_negatives': None,
    'precision': 0.576271186440678,
    'recall': 0.9927007299270073,
    'f1': 0.7292225201072386,
    'accuracy': 0.9966317614886947,
    'negative_predictive_value': 0.9999663865546219}