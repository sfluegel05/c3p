"""
Classifies: CHEBI:17761 ceramide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ceramide(smiles: str):
    """
    Determines if a molecule is a ceramide based on its SMILES string.
    Ceramides are characterized by an amide linkage with fatty acids, a sphingoid base, and typically have hydroxyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ceramide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"

    # Look for amide linkage pattern
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide linkage found"

    # Check for long fatty acid chain, typically C14-C26
    carbon_chain_pattern = Chem.MolFromSmarts("[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]")
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "No long carbon chain found"

    # Check for sphingoid base, nitrogen must be linked in sphingoid region
    sphingoid_base_pattern = Chem.MolFromSmarts("NC[C@H](O)C")
    if not mol.HasSubstructMatch(sphingoid_base_pattern):
        return False, "No sphingoid base structure found"

    # Check for hydroxyl groups at positions consistent with ceramides
    hydroxyl_pattern = Chem.MolFromSmarts("[C@H](O)")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 1:
        return False, "Missing hydroxyl group, which is common in ceramides"

    # Ensure overall structure supports ceramide classification
    n_atoms = mol.GetNumAtoms()
    n_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    n_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_nitrogens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)

    if n_nitrogens < 1 or n_carbons > n_atoms - n_oxygens - n_nitrogens:
        return False, "Structure inconsistent with ceramide"

    return True, "Structure contains features consistent with ceramide classification"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:17761',
                          'name': 'ceramide',
                          'definition': 'Ceramides (N-acyl-sphingoid bases) '
                                        'are a major subclass of sphingoid '
                                        'base derivatives with an amide-linked '
                                        'fatty acid. The fatty acids are '
                                        'typically saturated or '
                                        'monounsaturated with chain lengths '
                                        'from 14 to 26 carbon atoms; the '
                                        'presence of a hydroxyl group on '
                                        'carbon 2 is fairly common. Ceramides '
                                        'are generally precursors of more '
                                        'complex sphingolipids. In the '
                                        'illustrated generalised structure, '
                                        'R(1) = OH, OX (where X = acyl, '
                                        'glycosyl, phosphate, phosphonate, '
                                        'etc.), or H.',
                          'parents': ['CHEBI:26739', 'CHEBI:37622'],
                          'xrefs': [   'KEGG:C00195',
                                       'LIPID_MAPS_class:LMSP02',
                                       'PMID:7542630',
                                       'Wikipedia:Ceramide'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'C[C@@H]1OC(=O)C=C1',
                                     'name': '(S)-5-methylfuran-2(5H)-one',
                                     'reason': 'No amide linkage found'},
                                 {   'smiles': 'Cc1cccc(c1)C(O)=O',
                                     'name': 'm-toluic acid',
                                     'reason': 'No amide linkage found'},
                                 {   'smiles': 'CN(C)CC(=O)NC[C@@H]1[C@@H]([C@@H](N1)CO)C2=CC=C(C=C2)C3=CC=CC(=C3)C#N',
                                     'name': 'N-[[(2S,3S,4R)-3-[4-(3-cyanophenyl)phenyl]-4-(hydroxymethyl)-2-azetidinyl]methyl]-2-(dimethylamino)acetamide',
                                     'reason': 'No long carbon chain found'},
                                 {   'smiles': 'C[C@@H]1O[C@@H](OCCCCCCCCCC(=O)CC(O)=O)[C@H](O)C[C@H]1O',
                                     'name': 'bkos#20',
                                     'reason': 'No amide linkage found'},
                                 {   'smiles': 'C([C@@](COC(CCCCCCCCC/C=C\\CCCCCC)=O)(OC(CCCCCCCCC/C=C\\CCCCCC)=O)[H])OP(=O)(O)OC[C@@](CO)([H])O',
                                     'name': 'PG(18:1(11Z)/18:1(11Z))',
                                     'reason': 'No amide linkage found'},
                                 {   'smiles': 'O=C(N1C(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(CO)CC=2C3=C(C=CC=C3)NC2)CCC(=O)N)(C)C)(C)C)(CC)C)CCC1)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C)CC4=CC=CC=C4)(C)C)CO)(C)C)(C)C)CC(C)C)CCC(=O)N)(C)C)(C)C)C)C)(C)C',
                                     'name': 'Chrysospermin B',
                                     'reason': 'No long carbon chain found'},
                                 {   'smiles': 'CC1(C=CC2=C3C(=CC(=C2O1)OC)C(=O)C(=CO3)C4=CC5=C(C=C4OC)OCO5)C',
                                     'name': '6-methoxy-3-(6-methoxy-1,3-benzodioxol-5-yl)-8,8-dimethyl-4-pyrano[2,3-h][1]benzopyranone',
                                     'reason': 'No amide linkage found'},
                                 {   'smiles': 'O[C@@H]([C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)NCC(O)=O)C',
                                     'name': 'Glu-Thr-Gly',
                                     'reason': 'No long carbon chain found'},
                                 {   'smiles': 'O([C@H]1[C@H](O[C@@H]2O[C@H]([C@@H](O)[C@@H](O)[C@@H]2O)C)[C@@H](NC(=O)C)[C@@H](O[C@@H]1CO)O[C@H]3[C@@H](O)[C@H](O)[C@H](O[C@@H]3OC[C@H]4O[C@@H](O[C@H]5[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]5CO)O)[C@@H](O)[C@@H](O)[C@@H]4O)CO)[C@@H]6O[C@@H]([C@H](O)[C@H](O[C@]7(O[C@H]([C@H](NC(=O)C)[C@@H](O)C7)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]6O)CO',
                                     'name': '(2S,4S,5R,6R)-5-Acetamido-2-[(2S,3R,4S,5S,6R)-2-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2S,3S,4S,5S,6R)-2-[[(2R,3S,4S,5S,6S)-6-[(2R,3S,4R,5R)-5-acetamido-4,6-dihydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,4,5-trihydroxyoxan-2-yl]methoxy]-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]oxy-2-(hydroxymethyl)-4-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-3-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-4-yl]oxy-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'No long carbon chain found'},
                                 {   'smiles': 'COCCOC(C)C(O)=O',
                                     'name': '2-(2-methoxyethoxy)propanoic '
                                             'acid',
                                     'reason': 'No amide linkage found'}],
    'sample_false_negatives': [   {   'smiles': '[C@H]1([C@H](O[C@@H](OC[C@@H]([C@@H](CCCCCCCCCCCCCCC)O)OC(CCCCCCCCCCCCC/C=C\\CCCCCCCC)=O)[C@@H]([C@H]1O)O)CO)O',
                                      'name': "beta-D-glucosyl-(1<->1')-N-[(15Z)-tetracosenoyl]sphinganine",
                                      'reason': 'No amide linkage found'}],
    'attempt': 0,
    'reasoning': 'To classify chemical entities as ceramides based on their '
                 'SMILES string, we need to follow these steps:\n'
                 '\n'
                 '1. **Understand the Structure of Ceramides**: \n'
                 '   - Ceramides are primarily defined as having an '
                 'amide-linked fatty acid.\n'
                 '   - The fatty acids are typically saturated or '
                 'monounsaturated, with chain lengths from 14 to 26 carbon '
                 'atoms.\n'
                 '   - The presence of a hydroxyl group on carbon 2 is '
                 'common.\n'
                 '   - Ceramides can be precursors to more complex '
                 'sphingolipids.\n'
                 '\n'
                 '2. **SMILES Parsing and Pattern Identification**:\n'
                 '   - We will parse the SMILES string to identify key '
                 'structural features: the amide bond, long carbon chains '
                 'typical of fatty acids, and the presence of the sphingoid '
                 'base.\n'
                 '   - Look for hydroxyl groups and other structural features '
                 'consistent with ceramides.\n'
                 '\n'
                 '3. **Using RDKit for Structural Analysis**:\n'
                 '   - Utilize RDKit to search for these substructures and '
                 'calculate properties such as chain length.\n'
                 '\n'
                 '4. **Function Definition**:\n'
                 '   - Define a function `is_ceramide` that returns a boolean '
                 'value and a string explanation of the classification '
                 'result.\n'
                 '\n'
                 'Here is the Python code implementing this:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 570,
    'num_false_positives': 100,
    'num_true_negatives': 7330,
    'num_false_negatives': 1,
    'num_negatives': None,
    'precision': 0.8507462686567164,
    'recall': 0.9982486865148862,
    'f1': 0.9186140209508461,
    'accuracy': 0.987376577927759,
    'negative_predictive_value': 0.9998635929613968}