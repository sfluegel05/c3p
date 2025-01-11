"""
Classifies: CHEBI:36976 nucleotide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleotide(smiles: str):
    """
    Identifies if a molecule is a nucleotide based on its SMILES string.
    A nucleotide is a nucleoside phosphate, consisting of a nitrogenous base, sugar, and phosphate group.

    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is a nucleotide, False otherwise.
        str: Reason for classification.
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for different phosphate group occurrences
    phosphate_patterns = [
        Chem.MolFromSmarts('[OX1]P(=O)([OX1])[O-]'),  # Linear phosphates, common
        Chem.MolFromSmarts('[OX2]P(=O)([OX2])[OX2]'),  # Terminal phosphates
        Chem.MolFromSmarts('P([OX2H1])(=O)O[OX2]'),   # Cyclic phosphate
    ]
    
    phosphate_found = any(mol.HasSubstructMatch(pattern) for pattern in phosphate_patterns)
    if not phosphate_found:
        return False, "No phosphate group found"
    
    # Check for sugar rings (ribose and deoxyribose)
    ribose_patterns = [
        Chem.MolFromSmarts('C1CO[C@H](O)[C@@H]1O'),  # D-Ribose stereo and connectivity
        Chem.MolFromSmarts('C1COC(C1O)O'),           # Simplified ribose for catch-all
    ]
    deoxyribose_patterns = [
        Chem.MolFromSmarts('C1COC[C@H]1O'),          # D-Deoxyribose with connectivity
    ]

    sugar_found = any(mol.HasSubstructMatch(pattern) for pattern in ribose_patterns + deoxyribose_patterns)
    if not sugar_found:
        return False, "No compatible sugar ring found (ribose or deoxyribose)"
    
    # Check for nitrogenous bases
    base_patterns = [
        Chem.MolFromSmarts('n1cnc2c1ncnc2'),  # Purine base (common to adenine, guanine)
        Chem.MolFromSmarts('c1c[nH]c(=O)[nH]c1=O'),  # Pyrimidine base (uracil, thymine)
        Chem.MolFromSmarts('c1[nH]c(=O)[nH]c2c1ncnc2'),  # Cytosine variation
    ]

    base_found = any(mol.HasSubstructMatch(pattern) for pattern in base_patterns)
    if not base_found:
        return False, "No nitrogenous base found in the structure"

    return True, "Molecule matches nucleotide structure: base, sugar, and phosphate detected"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36976',
                          'name': 'nucleotide',
                          'definition': 'A nucleotide is a nucleoside '
                                        'phosphate resulting from the '
                                        'condensation of the 3 or 5 hydroxy '
                                        'group of a nucleoside with phosphoric '
                                        'acid.',
                          'parents': ['CHEBI:25608'],
                          'xrefs': ['KEGG:C00215', 'Wikipedia:Nucleotide'],
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
               'P(OC[C@H]1O[C@@H](N2C3=NC=NC(N)=C3N=C2)[C@H](O)[C@@H]1O)(OP(OC[C@H]4O[C@H]([N+]=5C=C(C=CC5)C(=O)N)[C@H](O)[C@@H]4O)([O-])=O)(O)=O '
               'NAME: alpha-Nicotinamide Adenine Dinucleotide REASON: MISSED '
               'No compatible sugar ring (ribose or deoxyribose) found\n'
               ' * SMILES: '
               'Nc1nc2n(cnc2c(=O)[nH]1)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OP(O)(O)=O)[C@@H](O)[C@H]1O '
               'NAME: GTP REASON: MISSED No phosphate group found\n'
               ' * SMILES: '
               'Nc1nc(=O)n(cc1CO)[C@H]1C[C@H](O)[C@@H](COP(O)(O)=O)O1 NAME: '
               '5-hydroxymethyldeoxycytidylic acid REASON: MISSED No phosphate '
               'group found\n'
               ' * SMILES: '
               'C[C@@H](OP(O)(=O)OC[C@@H](O)[C@@H](O)[C@@H](O)Cn1c2cc(C)c(C)cc2nc2c1nc(=O)[nH]c2=O)[C@H](N)C(O)=O '
               'NAME: FMN-L-threonine REASON: MISSED No phosphate group found\n'
               ' * SMILES: '
               'C1=C(CO[C@H]2[C@@H]([C@H]([C@@H]([C@H](O2)CO)O)O)O)C(NC(N1[C@@H]3O[C@H](COP(=O)(O)O)[C@H](C3)O)=O)=O '
               'NAME: (beta-D-glucopyranosyloxymethyl)deoxyuridine '
               "5'-monophosphate REASON: MISSED No phosphate group found\n"
               ' * SMILES: '
               'O[C@@H]1[C@@H](COP(O)(=O)OP(O)(=O)OP(O)(O)=O)O[C@H]([C@@H]1O)n1ccc(=O)[nH]c1=O '
               'NAME: UTP REASON: MISSED No phosphate group found\n'
               ' * SMILES: '
               'P(O[C@@H]1[C@H](O[C@H]([C@@H]1O)N2C=NC=3C(=NC=NC32)N)COP(O[C@@H]4[C@H](O[C@H]([C@@H]4O)N5C=NC=6C(NC(=NC65)N)=O)CO)(=O)[O-])([O-])(=O)[O-] '
               "NAME: G[3'-5']pAp[3'] REASON: MISSED No compatible sugar ring "
               '(ribose or deoxyribose) found\n'
               ' * SMILES: '
               'N1=C(C)C(NC(N1[C@@H]2O[C@H](COP(O)(=O)O)[C@H](C2)O)=O)=O NAME: '
               "6-azathymidine 5'-monophosphate REASON: MISSED No phosphate "
               'group found\n'
               ' * SMILES: '
               'O[C@@H]1[C@@H](COP(O)(=O)OP(O)(=O)OP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2ccc(=O)[nH]c2=O)O[C@H]([C@@H]1O)n1ccc(=O)[nH]c1=O '
               "NAME: P(1),P(4)-bis(uridin-5'-yl) tetraphosphate REASON: "
               'MISSED No phosphate group found\n'
               ' * SMILES: P1(OC2C(O1)C(OC2N3C4=NC=NC(N)=C4N=C3)CO)(O)=O NAME: '
               "Adenosine 2',3'-cyclic phosphate REASON: MISSED No phosphate "
               'group found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [],
    'sample_false_negatives': [   {   'smiles': 'Nc1nc(=O)n(cc1CO)[C@H]1C[C@H](O)[C@@H](COP(O)(O)=O)O1',
                                      'name': '5-hydroxymethyldeoxycytidylic '
                                              'acid',
                                      'reason': 'No nitrogenous base found in '
                                                'the structure'},
                                  {   'smiles': 'C[C@@H](OP(O)(=O)OC[C@@H](O)[C@@H](O)[C@@H](O)Cn1c2cc(C)c(C)cc2nc2c1nc(=O)[nH]c2=O)[C@H](N)C(O)=O',
                                      'name': 'FMN-L-threonine',
                                      'reason': 'No compatible sugar ring '
                                                'found (ribose or '
                                                'deoxyribose)'},
                                  {   'smiles': 'C1=C(CO[C@H]2[C@@H]([C@H]([C@@H]([C@H](O2)CO)O)O)O)C(NC(N1[C@@H]3O[C@H](COP(=O)(O)O)[C@H](C3)O)=O)=O',
                                      'name': '(beta-D-glucopyranosyloxymethyl)deoxyuridine '
                                              "5'-monophosphate",
                                      'reason': 'No nitrogenous base found in '
                                                'the structure'},
                                  {   'smiles': 'O[C@@H]1[C@@H](COP(O)(=O)OP(O)(=O)OP(O)(O)=O)O[C@H]([C@@H]1O)n1ccc(=O)[nH]c1=O',
                                      'name': 'UTP',
                                      'reason': 'No nitrogenous base found in '
                                                'the structure'},
                                  {   'smiles': 'P(O[C@@H]1[C@H](O[C@H]([C@@H]1O)N2C=NC=3C(=NC=NC32)N)COP(O[C@@H]4[C@H](O[C@H]([C@@H]4O)N5C=NC=6C(NC(=NC65)N)=O)CO)(=O)[O-])([O-])(=O)[O-]',
                                      'name': "G[3'-5']pAp[3']",
                                      'reason': 'No phosphate group found'},
                                  {   'smiles': 'N1=C(C)C(NC(N1[C@@H]2O[C@H](COP(O)(=O)O)[C@H](C2)O)=O)=O',
                                      'name': "6-azathymidine 5'-monophosphate",
                                      'reason': 'No nitrogenous base found in '
                                                'the structure'},
                                  {   'smiles': 'O[C@@H]1[C@@H](COP(O)(=O)OP(O)(=O)OP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2ccc(=O)[nH]c2=O)O[C@H]([C@@H]1O)n1ccc(=O)[nH]c1=O',
                                      'name': "P(1),P(4)-bis(uridin-5'-yl) "
                                              'tetraphosphate',
                                      'reason': 'No nitrogenous base found in '
                                                'the structure'},
                                  {   'smiles': 'O(C[C@H]1O[C@H]([C@@H]([C@@H]1O)O)N2C(=O)N=C(N)C=C2)P(O[C@@H]3[C@H](O[C@H]([C@@H]3O)N4C(=O)N=C(N)C=C4)COP(=O)(O)O)(=O)O',
                                      'name': 'pCpC',
                                      'reason': 'No nitrogenous base found in '
                                                'the structure'},
                                  {   'smiles': 'P1(O[C@]2([C@@](O1)(C(O[C@@H]2CO)N3C=CC(=NC3=O)N)[H])[H])(O)=O',
                                      'name': "Cytidine 2',3'-cyclic phosphate",
                                      'reason': 'No nitrogenous base found in '
                                                'the structure'},
                                  {   'smiles': 'C12=NC(NC(C1=NC=3C(N2C[C@@H]([C@@H]([C@@H](COP(=O)(O)O)O)O)O)=CC(=C(C3)C)C)=O)=O',
                                      'name': 'FMN',
                                      'reason': 'No compatible sugar ring '
                                                'found (ribose or '
                                                'deoxyribose)'}],
    'attempt': 3,
    'reasoning': 'Analyzing the previous program and the outcomes, the '
                 'function `is_nucleotide` attempted to identify nucleotides '
                 'based on the expected presence of phosphate groups, sugar '
                 '(ribose or deoxyribose rings), and nitrogenous bases. '
                 'However, it appears the SMILES patterns and algorithm could '
                 'not identify these components correctly, possibly due to '
                 'specificity or missing variations of these substructures. '
                 "Let's improve the approach:\n"
                 '\n'
                 '1. **Phosphate Group Detection:** The phosphate pattern '
                 'might be too narrow. Nucleotides can have various phosphate '
                 'forms (mono, di, tri, cyclic, linear) that need to be '
                 'detected.\n'
                 '\n'
                 '2. **Sugar Detection:** The SMARTS pattern used for '
                 'identifying ribose and deoxyribose sugars might lack '
                 'variations in possible stereochemistry or linkage patterns, '
                 'thus failing valid matches.\n'
                 '\n'
                 '3. **Nitrogenous Base Detection:** More comprehensive '
                 'patterns for nitrogenous bases, ensuring coverage for '
                 'structural variations.\n'
                 '\n'
                 '4. **Combining Verification:** After identifying parts '
                 'individually, we should ensure they are connected '
                 'appropriately to form a nucleotide: base-sugar-phosphate.\n'
                 '\n'
                 "With these considerations, let's refine the program:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 209,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 100,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 0.6763754045307443,
    'f1': 0.806949806949807,
    'accuracy': 0.6763754045307443,
    'negative_predictive_value': 0.0}