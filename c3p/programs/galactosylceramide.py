"""
Classifies: CHEBI:36498 galactosylceramide
"""
from rdkit import Chem

def is_galactosylceramide(smiles: str):
    """
    Determines if a molecule is a galactosylceramide based on its SMILES string.
    A galactosylceramide is a cerebroside with a galactose monosaccharide head group.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is a galactosylceramide, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Updated sphingosine or phytosphingosine backbone patterns
    sphingosine_patterns = [
        Chem.MolFromSmarts("N[C@@H](CO)COC1O[C@H](CO)[C@H](O)[C@@H]([C@@H]1O)O"),    # Sphingosine backbone
        Chem.MolFromSmarts("N[C@H](CO)COC1O[C@@H](CO)[C@H](O)[C@H](O)[C@@H]1O"),   # Phytosphere backbone
    ]
    
    if not any(mol.HasSubstructMatch(pattern) for pattern in sphingosine_patterns):
        return False, "No sphingosine backbone found"

    # Long hydrocarbon chains attached via amide linkage
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide linkage found"

    # Flexible pattern for galactose (including variations: with or without sulfate)
    galactose_patterns = [
        Chem.MolFromSmarts("CO[C@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@@H]1O"),         # beta-D-Galactose
        Chem.MolFromSmarts("CO[C@@H]1O[C@@H](CO)[C@H](O)[C@@H](O)[C@H]1O"),       # alpha-D-Galactose
        Chem.MolFromSmarts("C1[C@H](O)[C@@H](OS(=O)(=O)O)[C@H](O)[C@@H](CO)O1")   # Sulfo beta-D-Galactose
    ]
    
    if not any(mol.HasSubstructMatch(pattern) for pattern in galactose_patterns):
        return False, "No flexible galactose head group found"

    return True, "Contains a sphingosine backbone with long-chain amide linkage and a galactose head group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36498',
                          'name': 'galactosylceramide',
                          'definition': 'Any of the cerebrosides in which the '
                                        'monosaccharide head group is '
                                        'galactose.',
                          'parents': [   'CHEBI:23079',
                                         'CHEBI:5254',
                                         'CHEBI:62941'],
                          'xrefs': [   'KEGG:C02686',
                                       'KEGG:G11121',
                                       'PMID:16758576',
                                       'PMID:17855742',
                                       'PMID:2088646',
                                       'PMID:23065187',
                                       'PMID:23446636',
                                       'PMID:23650721',
                                       'PMID:25947378',
                                       'PMID:26058499',
                                       'PMID:26907993',
                                       'PMID:27786470'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: NONE\n'
               'False positives: SMILES: '
               'CCCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CCCCCCCCCCCCCCCC)CO[C@@H]1O[C@H](CO[C@@H]2O[C@@H](C)[C@@H](O)[C@@H](O)[C@@H]2O)[C@@H](O[C@@H]2O[C@H](CO)[C@@H](O[C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]3O)[C@H](O)[C@H]2NC(C)=O)[C@H](O)[C@H]1NC(C)=O '
               'NAME: CF4-C182L REASON: WRONGLY CLASSIFIED Contains a general '
               'sphingosine backbone with long-chain amide linkage and a '
               'flexible galactose head group\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCC(O)C(CO[C@H]1O[C@@H]([C@@H](O[C@H]2O[C@H](CO[C@H]3O[C@H](CO)[C@H](O)[C@H](O)[C@H]3O[C@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]3O)[C@@H](O)[C@H](O)[C@H]2N)[C@H](O)[C@H]1O)C(O)=O)NC(=O)[C@@H](O)CCCCCCCCCCCC '
               'NAME: '
               '3-hydroxy-2-{[(2S)-2-hydroxytetradecanoyl]amino}octadecyl '
               'alpha-D-Man-(1->2)-alpha-D-Gal-(1->6)-alpha-D-GlcN-(1->4)-alpha-D-GlcA '
               'REASON: WRONGLY CLASSIFIED Contains a general sphingosine '
               'backbone with long-chain amide linkage and a flexible '
               'galactose head group\n'
               ' * SMILES: '
               'CCCCCCCCCCCC[C@H](O)C(=O)NC(CO[C@H]1O[C@@H]([C@@H](O[C@H]2O[C@H](CO[C@H]3O[C@H](CO)[C@H](O)[C@H](O)[C@H]3O[C@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]3O)[C@@H](O)[C@H](O)[C@H]2N)[C@H](O)[C@H]1O)C(O)=O)C(O)CCCCCCCCCC1CC1CCCCCC '
               'NAME: '
               '12-(2-hexylcyclopropyl)-3-hydroxy-2-{[(2S)-2-hydroxytetradecanoyl]amino}dodecyl '
               'alpha-D-Man-(1->2)-alpha-D-Gal-(1->6)-alpha-D-GlcN-(1->4)-alpha-D-GlcA '
               'REASON: WRONGLY CLASSIFIED Contains a general sphingosine '
               'backbone with long-chain amide linkage and a flexible '
               'galactose head group\n'
               ' * SMILES: '
               'O([C@@H]1[C@H](O[C@@H](OC[C@@H]([C@@H](CCCCCCCCCCCCCCC)O)NC(CCCCCCCCCCCCCCCCC)=O)[C@@H]([C@H]1O)O)CO)[C@H]2[C@@H]([C@H]([C@H]([C@H](O2)CO)O[C@H]3[C@@H]([C@H]([C@H]([C@H](O3)CO)O)O[C@H]4[C@@H]([C@H]([C@H]([C@H](O4)CO)O)O)O)NC(C)=O)O[C@@]5(C[C@@H]([C@H]([C@@](O5)([C@@H]([C@@H](CO)O)O)[H])NC(C)=O)O)C(=O)O)O '
               'NAME: '
               "beta-D-Gal-(1->3)-beta-D-GalNAc-(1->4)-[alpha-Neu5Ac-(2->3)]-beta-D-Gal-(1->4)-beta-D-Glc-(1<->1')-Cer(18:0/18:0) "
               'REASON: WRONGLY CLASSIFIED Contains a general sphingosine '
               'backbone with long-chain amide linkage and a flexible '
               'galactose head group\n'
               'False negatives: SMILES: '
               'C(=C/CCCCCCCCCCCCC)\\[C@@H](O)[C@@H](NC(=O)CCCCCCC/C=C\\CCCCCCCC)COC1[C@@H]([C@H]([C@H]([C@H](O1)CO)O)O)O '
               'NAME: 1-(beta-D-galactosyl)-N-(9Z-octadecenoyl)-sphingosine '
               'REASON: MISSED No general sphingosine backbone found\n'
               ' * SMILES: '
               '[C@H]1([C@@H]([C@H]([C@H]([C@H](O1)CO)O)O)O)OC[C@@H]([C@@H]([C@@H](CCCCCCCC=2C=CC=CC2)O)O)NC(CCCCCCCCCCCCCCCCCCCCCCCCC)=O '
               'NAME: '
               'N-[(2S,3S,4R)-1-(alpha-D-galactosyloxy)-3,4-dihydroxy-11-phenylundecan-2-yl]hexacosanamide '
               'REASON: MISSED No general sphingosine backbone found\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](OS(O)(=O)=O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC '
               'NAME: 1-(3-O-sulfo-beta-D-galactosyl)-N-palmitoylsphingosine '
               'REASON: MISSED No general sphingosine backbone found\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCCCC[C@@H](O)C(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC '
               'NAME: '
               '1-(beta-D-galactosyl)-N-[(2R)-2-hydroxybehenoyl]sphingosine '
               'REASON: MISSED No general sphingosine backbone found\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](COC1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC '
               'NAME: N-lignoceroyl-D-galactosylsphingosine REASON: MISSED No '
               'general sphingosine backbone found\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC '
               'NAME: beta-D-galactosyl-N-eicosanoylsphingosine REASON: MISSED '
               'No general sphingosine backbone found\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCCCCCCC(O)C(=O)N[C@@H](COC1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC '
               'NAME: N-(2-hydroxylignoceroyl)-D-galactosylsphingosine REASON: '
               'MISSED No general sphingosine backbone found\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC '
               'NAME: beta-D-galactosyl-N-octadecanoylsphingosine REASON: '
               'MISSED No general sphingosine backbone found\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](COC1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC '
               'NAME: N-tricosanoyl-D-galactosylsphingosine REASON: MISSED No '
               'general sphingosine backbone found\n'
               ' * SMILES: '
               '[C@@H]1([C@H](O[C@@H](OC[C@@H]([C@@H](/C=C/CCCCCCCCCCCCC)O)NC(CCCCCCCCCCCCC/C=C\\CCCCCCCC)=O)[C@@H]([C@H]1O)O)CO)O '
               'NAME: beta-D-galactosyl-N-[(15Z)-tetracosenoyl]sphingosine '
               'REASON: MISSED No general sphingosine backbone found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'O([C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)CO)[C@H]3[C@@H](O)[C@H](OC(O)[C@@H]3O)CO',
                                     'name': 'beta-D-Glcp-(1->2)-beta-D-Glcp-(1->3)-D-Galp',
                                     'reason': 'No sphingosine backbone found'},
                                 {   'smiles': 'CCS(=O)(=O)NCC[C@@H]1CC[C@@H]([C@H](O1)CO)NC(=O)NC2=CC(=CC(=C2)Cl)Cl',
                                     'name': '1-(3,5-dichlorophenyl)-3-[(2S,3S,6S)-6-[2-(ethylsulfonylamino)ethyl]-2-(hydroxymethyl)-3-oxanyl]urea',
                                     'reason': 'No sphingosine backbone found'},
                                 {   'smiles': 'C(C(O)=O)C/C=C\\C/C=C\\C\\C=C/C=C/C=C/[C@H]1[C@H](C/C=C\\CC)O1',
                                     'name': '(16S,17S)-epoxy-(4Z,7Z,10Z,12E,14E,19Z)-docosahexaenoic '
                                             'acid',
                                     'reason': 'No sphingosine backbone found'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@H]2[C@@H](CO)O[C@@H](O[C@@H]3[C@@H](CO)O[C@@H](O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O',
                                     'name': 'alpha-D-Galp-(1->4)-beta-D-Galp-(1->4)-beta-D-Glcp',
                                     'reason': 'No sphingosine backbone found'},
                                 {   'smiles': '[Li+].[Br-]',
                                     'name': 'lithium bromide',
                                     'reason': 'No sphingosine backbone found'},
                                 {   'smiles': 'C=1(OC)C2=C(C=C3C[C@H]([C@](CC=4C=C(OC)C(OC)=C(C4C13)OC)(C)O)C)OCO2',
                                     'name': 'Besigomsin',
                                     'reason': 'No sphingosine backbone found'},
                                 {   'smiles': 'C1=CC=C(C(=C1)C=CC(=O)C2=CC=CN2)Cl',
                                     'name': '3-(2-chlorophenyl)-1-(1H-pyrrol-2-yl)-2-propen-1-one',
                                     'reason': 'No sphingosine backbone found'},
                                 {   'smiles': 'O=C1OC(O)C(=C1C(O)C(C)C)C',
                                     'name': '5-Hydroxy-3-(1-hydroxy-2-methylpropyl)-4-methyl-2(5H)-furanone',
                                     'reason': 'No sphingosine backbone found'},
                                 {   'smiles': '[H][C@@]1(COC(C)(C)[C@@]1([H])OC(=O)\\C=C\\c1ccccc1)c1c(O)ccc2C(=O)C[C@H](Oc12)c1ccccc1',
                                     'name': '(+)-tephrorin B',
                                     'reason': 'No sphingosine backbone found'},
                                 {   'smiles': 'C[C@@H]1CN(C(=O)C2=C(C3=CC=CC=C3CO[C@@H]1CN(C)CC4=CC=CC(=C4)C(=O)O)C5=CC=CC=C5N2C)[C@@H](C)CO',
                                     'name': 'LSM-9341',
                                     'reason': 'No sphingosine backbone '
                                               'found'}],
    'sample_false_negatives': [],
    'attempt': 5,
    'reasoning': 'To improve the classification of galactosylceramide based on '
                 'their SMILES strings, we need to address the shortcomings of '
                 "the previous implementation. Here's a step-by-step "
                 'hypothesis and plan for improvement:\n'
                 '\n'
                 '1. **Identify Critical Patterns**: \n'
                 '   - We need to ensure we capture the presence of the '
                 'specific ceramide components and recognize variations like '
                 'hydroxylation and sulfate groups.\n'
                 '   - The previous code was too narrow in terms of defining '
                 'the galactose head group and sphingosine backbone.\n'
                 '\n'
                 '2. **Adjust Substructure Patterns**:\n'
                 '   - Improve the SMARTS patterns to better capture the sugar '
                 'moiety. Galactosylceramides may have alpha or beta linkages '
                 'and possible sulfations.\n'
                 '   - Ensure that the pattern correctly identifies '
                 'sphingosine or phytosphingosine backbones, considering '
                 'common hydroxylations and unsaturations.\n'
                 '\n'
                 '3. **Use Multiple Patterns**:\n'
                 '   - Incorporate multiple substructure searches to ensure '
                 'that alternative structures and chiral variations are '
                 'considered.\n'
                 '\n'
                 '4. **Check Results**:\n'
                 '   - Verify and test SMILES strings known to correspond to '
                 'the galactosylceramide class.\n'
                 '\n'
                 "Let's update the function accordingly:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 37,
    'num_false_positives': 100,
    'num_true_negatives': 5101,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.27007299270072993,
    'recall': 1.0,
    'f1': 0.42528735632183906,
    'accuracy': 0.9809087437953418,
    'negative_predictive_value': 1.0}