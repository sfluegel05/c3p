"""
Classifies: CHEBI:47622 acetate ester
"""
from rdkit import Chem

def is_acetate_ester(smiles: str):
    """
    Determines if a molecule is an acetate ester based on its SMILES string.
    An acetate ester contains the structural pattern of CH3-C(=O)-O-.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an acetate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define acetate ester pattern (CH3-C(=O)-O-)
    acetate_pattern = Chem.MolFromSmarts("CC(=O)O")
    if acetate_pattern is None:
        return False, "Error in SMARTS pattern creation"

    # Search for the acetate pattern in the molecule
    if mol.HasSubstructMatch(acetate_pattern):
        return True, "Contains acetate ester functional group"
    else:
        return False, "Does not contain acetate ester structural pattern"

# Testing some examples
test_smiles = [
    "CC(=O)OC",    # ethyl acetate
    "CCCCCOC(C)=O",    # pentyl acetate
    "O=C(C)OC",    # methyl acetate
    "CC1=CC(=O)OC1",    # acetoxyethane
]

for smiles in test_smiles:
    result, reason = is_acetate_ester(smiles)
    print(f"SMILES: {smiles} -> Acetate Ester: {result}, Reason: {reason}")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:47622',
                          'name': 'acetate ester',
                          'definition': 'Any carboxylic ester where the '
                                        'carboxylic acid component is acetic '
                                        'acid.',
                          'parents': ['CHEBI:33308'],
                          'xrefs': ['KEGG:C01883', 'Wikipedia:Acetate#Esters'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'O=C1C2=C(OC(=C1)C)C3=C(OC)C=C(OC)C=C3C(=C2O)C4=C5OC(=CC(C5=C(O)C=6C4=CC(OC)=CC6OC)=O)C',
                                     'name': 'Isonigerone',
                                     'reason': 'Does not contain acetate ester '
                                               'structural pattern'},
                                 {   'smiles': 'C(=C\\C/C=C\\CCCC(NC)=O)\\C/C=C\\C/C=C\\CCCCC',
                                     'name': 'N-methyl arachidonoyl amine',
                                     'reason': 'Does not contain acetate ester '
                                               'structural pattern'},
                                 {   'smiles': 'OC1(O)[C@]23N(CC1)C(N(O)[C@H]([C@@]3(N=C(N2)N)[H])CO)=N',
                                     'name': 'Decarbamoylneosaxitoxin',
                                     'reason': 'Does not contain acetate ester '
                                               'structural pattern'},
                                 {   'smiles': 'S([C@H]1N(C(=O)/C(=C\\C2=CC=C(OCC=C(C)C)C=C2)/N(C1=O)C)C)C',
                                     'name': 'Fusaperazine F',
                                     'reason': 'Does not contain acetate ester '
                                               'structural pattern'},
                                 {   'smiles': 'Oc1ccccc1I',
                                     'name': '2-iodophenol',
                                     'reason': 'Does not contain acetate ester '
                                               'structural pattern'},
                                 {   'smiles': 'O(C1=C(OC)C=C(C2=C(O)C(OC)=C(C3=CC=CC=C3)C=C2OC)C=C1)C/C=C(/CO)\\C',
                                     'name': 'Prenylterphenyllin F',
                                     'reason': 'Does not contain acetate ester '
                                               'structural pattern'},
                                 {   'smiles': 'O1[C@@H]([C@@H](O)[C@H](O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)[C@@H](O)[C@@H]1OC[C@H]3O[C@@H](OC[C@H]4O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]4O)[C@H](O)[C@@H](O)[C@@H]3O)CO[C@@H]5O[C@@H]([C@@H](O)[C@H](O[C@@H]6O[C@@H]([C@@H](O)[C@H](O)[C@H]6O)CO)[C@H]5O)CO[C@@H]7O[C@@H]([C@@H](O)[C@H](O)[C@H]7O)CO',
                                     'name': '(2R,3R,4S,5S,6R)-6-[[(2R,3R,4S,5S,6R)-6-[[(2R,3R,4S,5R,6R)-6-[[(2R,3R,4S,5R,6R)-3,5-Dihydroxy-4-[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-6-[[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]oxan-2-yl]oxymethyl]-3,5-dihydroxy-4-[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxymethyl]-3,4,5-trihydroxyoxan-2-yl]oxymethyl]oxane-2,3,4,5-tetrol',
                                     'reason': 'Does not contain acetate ester '
                                               'structural pattern'},
                                 {   'smiles': 'O([C@H]1[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]1CO)O[C@H]2[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]2CO[C@@H]3O[C@H]([C@@H](O)[C@@H](O)[C@@H]3O)C)O)[C@@H]4O[C@@H]([C@@H](O)[C@H](O[C@H]5O[C@@H]([C@@H](O)[C@H](O)[C@@H]5O[C@@H]6O[C@@H]([C@@H](O[C@@H]7O[C@@H]([C@H](O)[C@H](O)[C@H]7O)CO)[C@H](O)[C@@H]6NC(=O)C)CO)CO)[C@@H]4O)CO[C@H]8O[C@@H]([C@@H](O)[C@H](O)[C@@H]8O[C@@H]9O[C@@H]([C@@H](O[C@@H]%10O[C@@H]([C@H](O)[C@H](O)[C@H]%10O)CO)[C@H](O)[C@H]9NC(=O)C)CO)CO',
                                     'name': 'Gal2GlcNAc2Man3GlcNAcFucGlcNAc',
                                     'reason': 'Does not contain acetate ester '
                                               'structural pattern'},
                                 {   'smiles': 'CN1CC2(CCN(CC2)S(=O)(=O)C3=CN(C=N3)C)C4=C([C@@H]1CO)NC5=C4C=CC(=C5)OC',
                                     'name': "[(1R)-7-methoxy-2-methyl-1'-[(1-methyl-4-imidazolyl)sulfonyl]-1-spiro[3,9-dihydro-1H-pyrido[3,4-b]indole-4,4'-piperidine]yl]methanol",
                                     'reason': 'Does not contain acetate ester '
                                               'structural pattern'},
                                 {   'smiles': 'O=C(N(CCO)CCO)CN1C(=NC=C1)[N+]([O-])=O',
                                     'name': '1H-Imidazole-1-acetamide, '
                                             'N,N-bis(2-hydroxyethyl)-2-nitro-',
                                     'reason': 'Does not contain acetate ester '
                                               'structural pattern'}],
    'sample_false_negatives': [],
    'attempt': 0,
    'reasoning': 'To classify chemical entities as acetate esters based on '
                 'their SMILES strings, we need to identify the structural '
                 'features typical of an acetate ester. An acetate ester is '
                 'any carboxylic ester derived from acetic acid, which has the '
                 'substructural pattern of an acyl group (C(=O)O) bonded to an '
                 'acedic acid moiety. The structural hallmark of this class '
                 'would be the presence of an acetic acid ester linkage '
                 '(CH3-C(=O)-O-).\n'
                 '\n'
                 'Here is a step-by-step approach:\n'
                 '\n'
                 '1. **Import necessary modules from rdkit**: We will require '
                 'RDKit modules to parse and analyze the SMILES string.\n'
                 '\n'
                 '2. **Function Definition**:\n'
                 '   - The function `is_acetate_ester` will accept a SMILES '
                 'string as input.\n'
                 '   - It will return whether or not the string represents an '
                 'acetate ester, along with a reason for its conclusion.\n'
                 '\n'
                 '3. **Parse SMILES String**:\n'
                 '   - Use RDKit to convert the SMILES string into a molecular '
                 'object.\n'
                 '\n'
                 '4. **Define SMARTS Pattern**:\n'
                 '   - Use a SMARTS pattern to recognize the acetic acid ester '
                 'group typically represented as `CC(=O)O`.\n'
                 '\n'
                 '5. **Perform Matching**:\n'
                 '   - Check the molecular structure for the presence of the '
                 'acetic acid ester pattern.\n'
                 '   - Return the classification result with a reason.\n'
                 '\n'
                 "Now, let's write the code for the above logic:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 554,
    'num_false_positives': 100,
    'num_true_negatives': 176,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.8470948012232415,
    'recall': 1.0,
    'f1': 0.9172185430463576,
    'accuracy': 0.8795180722891566,
    'negative_predictive_value': 1.0}