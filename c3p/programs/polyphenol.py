"""
Classifies: CHEBI:26195 polyphenol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyphenol(smiles: str):
    """
    Determines if a molecule is a polyphenol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyphenol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for multiple aromatic rings containing hydroxyl groups
    phenol_pattern = Chem.MolFromSmarts("c1cc(O)ccc1")  # Phenol or aromatic hydroxyl group
    matches = mol.GetSubstructMatches(phenol_pattern)
    if len(matches) < 2:
        return False, f"Found {len(matches)} phenolic rings, need at least 2 for polyphenol"

    # Further confirm by counting the total number of aromatic rings
    num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if num_aromatic_rings < 2:
        return False, f"Found {num_aromatic_rings} aromatic rings, need at least 2 for polyphenol"

    return True, "Contains multiple phenolic rings, consistent with a polyphenol"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26195',
                          'name': 'polyphenol',
                          'definition': 'Members of the class of phenols that '
                                        'contain 2 or more benzene rings each '
                                        'of which is substituted by at least '
                                        'one hydroxy group.',
                          'parents': ['CHEBI:33853'],
                          'xrefs': ['Wikipedia:Polyphenol'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'P(OCC(OC(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCC/C=C\\C/C=C\\C/C=C\\CCCCCCCC)(OCCN(C)C)(O)=O',
                                     'name': 'Pe-nme2(20:3(5Z,8Z,11Z)/22:6(4Z,7Z,10Z,13Z,16Z,19Z))',
                                     'reason': 'Found 0 phenolic rings, need '
                                               'at least 2 for polyphenol'},
                                 {   'smiles': 'O1[C@]2([C@@H](O)[C@H](O)[C@@H](O)[C@@]1(OCC[C@H](CCC=C(C(OC[C@]3(O)[C@@H](O)[C@@](OC3)(OC2)[H])=O)C)C)[H])[H]',
                                     'name': 'Urceolide',
                                     'reason': 'Found 0 phenolic rings, need '
                                               'at least 2 for polyphenol'},
                                 {   'smiles': 'CC1=CC(=C(C=C1)C)C(=O)CSC2=NN=C(S2)C',
                                     'name': '1-(2,5-dimethylphenyl)-2-[(5-methyl-1,3,4-thiadiazol-2-yl)thio]ethanone',
                                     'reason': 'Found 0 phenolic rings, need '
                                               'at least 2 for polyphenol'},
                                 {   'smiles': 'O=C1C=C([C@@]2(C(C[C@@H](C2)O)(C)C)C)CC[C@@]1(O)C',
                                     'name': 'Enokipodin H',
                                     'reason': 'Found 0 phenolic rings, need '
                                               'at least 2 for polyphenol'},
                                 {   'smiles': 'CCCC[C@](Cn1cncn1)(C#N)c1ccc(Cl)cc1',
                                     'name': '(S)-myclobutanil',
                                     'reason': 'Found 0 phenolic rings, need '
                                               'at least 2 for polyphenol'},
                                 {   'smiles': 'O=C1O[C@@H]([C@H](NC(=O)[C@@H](NC(=O)CCCCC)CC(=O)O)C(=O)N[C@H](C(=O)N[C@H]2CC[C@H](N([C@H](C(N([C@H](C(N[C@H]1C(C)C)=O)CC3=CC=CC=C3)C)=O)CC(C)C)C2=O)O)CCCCN(C)C)C',
                                     'name': 'Cyanopeptolin D',
                                     'reason': 'Found 0 phenolic rings, need '
                                               'at least 2 for polyphenol'},
                                 {   'smiles': 'O1[C@]2(O)N([C@H](C(=O)N3[C@]2(CCC3)[H])C(C)C)C(=O)[C@@]1(NC(=O)[C@H]4CN([C@]5(C(=C4)C6=C7C(C5)=CNC7=CC=C6)[H])C)C',
                                     'name': 'Ergovaline',
                                     'reason': 'Found 0 phenolic rings, need '
                                               'at least 2 for polyphenol'},
                                 {   'smiles': 'OC(=O)CCCCCCC#CCCCC',
                                     'name': '8-tridecynoic acid',
                                     'reason': 'Found 0 phenolic rings, need '
                                               'at least 2 for polyphenol'},
                                 {   'smiles': 'CC(=O)N([O-])CCC[C@H](N)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CO)C(=O)N[C@H]([C@H](O)[C@H]1S[C@H]([C@H](O)[C@H]1O)n1ccc(=N)n(C)c1=O)C(O)=O',
                                     'name': 'desferrialbomycin epsilon(3-)',
                                     'reason': 'Found 0 phenolic rings, need '
                                               'at least 2 for polyphenol'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@H](O[C@@H]3[C@@H](CO)O[C@H](O[C@@H]4[C@@H](CO)O[C@H](O)[C@H](O)[C@H]4O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O',
                                     'name': 'alpha-maltotetraose',
                                     'reason': 'Found 0 phenolic rings, need '
                                               'at least 2 for polyphenol'}],
    'sample_false_negatives': [   {   'smiles': 'CC(=O)C1=C(O)C=C2Oc3c(c(O)cc(O)c3C(N)=O)[C@]2(C)C1=O',
                                      'name': 'cercosporamide',
                                      'reason': 'Found 1 aromatic rings, need '
                                                'at least 2 for polyphenol'},
                                  {   'smiles': 'CCCC(=O)c1c(O)c(C)c(O)c(CC2=C(O)C(C)(C)C(O)=C(C(C)=O)C2=O)c1O',
                                      'name': 'flavaspidic acid AB',
                                      'reason': 'Found 1 aromatic rings, need '
                                                'at least 2 for polyphenol'},
                                  {   'smiles': 'O=C(O)C1=CC(OC)=C(O[C@@H]2O[C@@H]([C@@H](OC)[C@@H]([C@H]2O)O)CO)C(=C1)CC=C(C)C',
                                      'name': 'Conoideoglucoside C',
                                      'reason': 'Found 1 aromatic rings, need '
                                                'at least 2 for polyphenol'},
                                  {   'smiles': 'COc1cc(ccc1O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)C(O)=O',
                                      'name': '4-(beta-D-glucopyranosyloxy)-3-methoxybenzoic '
                                              'acid',
                                      'reason': 'Found 1 aromatic rings, need '
                                                'at least 2 for polyphenol'},
                                  {   'smiles': 'C\\C=C\\c1c(O)c(O)cc2OC(=O)[C@](C)(NC(=O)\\C=C\\C(O)=O)c12',
                                      'name': 'fumimycin',
                                      'reason': 'Found 1 aromatic rings, need '
                                                'at least 2 for polyphenol'},
                                  {   'smiles': 'O1C(C(O)C(O)C(O)C1OC(=O)C=2C=C(OC)C(OC)=C(O)C2)CO',
                                      'name': '[3,4,5-Trihydroxy-6-(hydroxymethyl)oxan-2-yl] '
                                              '3-hydroxy-4,5-dimethoxybenzoate',
                                      'reason': 'Found 1 aromatic rings, need '
                                                'at least 2 for polyphenol'},
                                  {   'smiles': 'O1C(C(O)C(OC(=O)C=2C=C(OC3OC(C(O)C(O)C3O)C(O)=O)C(O)=C(O)C2)C(O)C1O)CO',
                                      'name': '6-[2,3-dihydroxy-5-({[2,3,5-trihydroxy-6-(hydroxymethyl)oxan-4-yl]oxy}carbonyl)phenoxy]-3,4,5-trihydroxyoxane-2-carboxylic '
                                              'acid',
                                      'reason': 'Found 1 aromatic rings, need '
                                                'at least 2 for polyphenol'},
                                  {   'smiles': 'O1C(C(O)C(O)C(O)C1OC=2C(O)=C(C=C(C2)C(OC)=O)CC=C(C)C)CO',
                                      'name': 'Methyl '
                                              '3,4-dihydroxy-5-prenylbenzoate '
                                              '3-glucoside',
                                      'reason': 'Found 1 aromatic rings, need '
                                                'at least 2 for polyphenol'},
                                  {   'smiles': 'O[C@@H]1C[C@](O)(C[C@@H](OC(=O)\\C=C\\c2ccc(O)c(O)c2)[C@@H]1O)C(O)=O',
                                      'name': 'chlorogenic acid',
                                      'reason': 'Found 1 aromatic rings, need '
                                                'at least 2 for polyphenol'},
                                  {   'smiles': 'O1C(C(O)C(O)C(O)C1OC=2C(O)=C(OC)C=C(C2)C(O)=O)C(O)=O',
                                      'name': '6-(5-carboxy-2-hydroxy-3-methoxyphenoxy)-3,4,5-trihydroxyoxane-2-carboxylic '
                                              'acid',
                                      'reason': 'Found 1 aromatic rings, need '
                                                'at least 2 for polyphenol'}],
    'attempt': 0,
    'reasoning': 'To create a program to classify polyphenols based on their '
                 'SMILES string, we need to identify the structural '
                 'characteristics that define a polyphenol. Polyphenols are '
                 'typically characterized by the presence of multiple phenol '
                 'units, often interconnected through various linkages, which '
                 'can include carbon-carbon bonds or other types of bonds. '
                 "Given the definition provided as 'None', we will rely on "
                 'common features of polyphenols for our classification.\n'
                 '\n'
                 'Polyphenol characteristics:\n'
                 '- Multiple aromatic (phenolic) rings.\n'
                 '- Hydroxyl groups attached to the aromatic rings.\n'
                 '- Often linked to form complex structures.\n'
                 '\n'
                 'For classification using SMILES, we need to identify '
                 'patterns that include these characteristics. We can use '
                 'SMARTS pattern matching to find the presence of these '
                 'features.\n'
                 '\n'
                 "Here's how we can approach this:\n"
                 '1. Parse the SMILES string to create a molecular structure.\n'
                 '2. Check for the presence of multiple aromatic rings with '
                 'hydroxyl groups.\n'
                 '3. Return true if these features are present with a suitable '
                 'explanation, otherwise return false.\n'
                 '\n'
                 "Let's implement this logic in Python using RDKit.",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 654,
    'num_false_positives': 100,
    'num_true_negatives': 961,
    'num_false_negatives': 28,
    'num_negatives': None,
    'precision': 0.8673740053050398,
    'recall': 0.9589442815249267,
    'f1': 0.9108635097493036,
    'accuracy': 0.9265633964429145,
    'negative_predictive_value': 0.9716885743174924}