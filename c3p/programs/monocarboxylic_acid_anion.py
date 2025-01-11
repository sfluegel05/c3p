"""
Classifies: CHEBI:35757 monocarboxylic acid anion
"""
from rdkit import Chem

def is_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a monocarboxylic acid anion based on its SMILES string.
    A monocarboxylic acid anion is formed when the carboxy group of a monocarboxylic acid is deprotonated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monocarboxylic acid anion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for carboxylate group [O-]C(=O)
    carboxylate_pattern = Chem.MolFromSmarts("[O-]C(=O)")
    
    # Find matches for carboxylate pattern
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    if len(carboxylate_matches) == 0:
        return False, "No carboxylate group found"
    if len(carboxylate_matches) > 1:
        return False, f"Multiple ({len(carboxylate_matches)}) carboxylate groups found"

    return True, "Contains a single deprotonated carboxylic acid group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35757',
                          'name': 'monocarboxylic acid anion',
                          'definition': 'A carboxylic acid anion formed when '
                                        'the carboxy group of a monocarboxylic '
                                        'acid is deprotonated.',
                          'parents': ['CHEBI:29067'],
                          'xrefs': ['KEGG:C00060'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'O([C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)CO)[C@H]3[C@@H](O)[C@H](OC(O)[C@@H]3O)CO',
                                     'name': 'beta-D-Glcp-(1->2)-beta-D-Glcp-(1->3)-D-Galp',
                                     'reason': 'No carboxylate group found'},
                                 {   'smiles': 'CCS(=O)(=O)NCC[C@@H]1CC[C@@H]([C@H](O1)CO)NC(=O)NC2=CC(=CC(=C2)Cl)Cl',
                                     'name': '1-(3,5-dichlorophenyl)-3-[(2S,3S,6S)-6-[2-(ethylsulfonylamino)ethyl]-2-(hydroxymethyl)-3-oxanyl]urea',
                                     'reason': 'No carboxylate group found'},
                                 {   'smiles': 'C(C(O)=O)C/C=C\\C/C=C\\C\\C=C/C=C/C=C/[C@H]1[C@H](C/C=C\\CC)O1',
                                     'name': '(16S,17S)-epoxy-(4Z,7Z,10Z,12E,14E,19Z)-docosahexaenoic '
                                             'acid',
                                     'reason': 'No carboxylate group found'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@H]2[C@@H](CO)O[C@@H](O[C@@H]3[C@@H](CO)O[C@@H](O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O',
                                     'name': 'alpha-D-Galp-(1->4)-beta-D-Galp-(1->4)-beta-D-Glcp',
                                     'reason': 'No carboxylate group found'},
                                 {   'smiles': '[Li+].[Br-]',
                                     'name': 'lithium bromide',
                                     'reason': 'No carboxylate group found'},
                                 {   'smiles': 'C=1(OC)C2=C(C=C3C[C@H]([C@](CC=4C=C(OC)C(OC)=C(C4C13)OC)(C)O)C)OCO2',
                                     'name': 'Besigomsin',
                                     'reason': 'No carboxylate group found'},
                                 {   'smiles': 'C1=CC=C(C(=C1)C=CC(=O)C2=CC=CN2)Cl',
                                     'name': '3-(2-chlorophenyl)-1-(1H-pyrrol-2-yl)-2-propen-1-one',
                                     'reason': 'No carboxylate group found'},
                                 {   'smiles': 'O=C1OC(O)C(=C1C(O)C(C)C)C',
                                     'name': '5-Hydroxy-3-(1-hydroxy-2-methylpropyl)-4-methyl-2(5H)-furanone',
                                     'reason': 'No carboxylate group found'},
                                 {   'smiles': '[H][C@@]1(COC(C)(C)[C@@]1([H])OC(=O)\\C=C\\c1ccccc1)c1c(O)ccc2C(=O)C[C@H](Oc12)c1ccccc1',
                                     'name': '(+)-tephrorin B',
                                     'reason': 'No carboxylate group found'},
                                 {   'smiles': 'C[C@@H]1CN(C(=O)C2=C(C3=CC=CC=C3CO[C@@H]1CN(C)CC4=CC=CC(=C4)C(=O)O)C5=CC=CC=C5N2C)[C@@H](C)CO',
                                     'name': 'LSM-9341',
                                     'reason': 'No carboxylate group found'}],
    'sample_false_negatives': [   {   'smiles': '[H][C@]12SC(C)(C)[C@@H](N1C(=O)[C@H]2NC(=O)C(C([O-])=O)c1ccccc1)C([O-])=O',
                                      'name': 'carbenicillin(2-)',
                                      'reason': 'Multiple (2) carboxylate '
                                                'groups found'},
                                  {   'smiles': 'S([O-])(=O)(=O)CCN(C(=O)CC[C@H]([C@@]1([C@@]2([C@]([C@]3([C@@]([C@@]4([C@]([C@@H](O)C3)(C[C@H](O)CC4)[H])C)(CC2)[H])[H])(CC1)[H])C)[H])C)[H]',
                                      'name': 'taurohyodeoxycholate(1-)',
                                      'reason': 'No carboxylate group found'},
                                  {   'smiles': 'C1=C[C@H](C=C[C@]1(C([O-])=O)CC(C(=O)[O-])=O)[NH3+]',
                                      'name': '4-amino-4-deoxyprephenate',
                                      'reason': 'Multiple (2) carboxylate '
                                                'groups found'},
                                  {   'smiles': 'CO[C@H]1O[C@H](COS([O-])(=O)=O)[C@@H](O[C@@H]2O[C@H]([C@@H](O[C@H]3O[C@H](COS([O-])(=O)=O)[C@@H](O[C@@H]4O[C@@H]([C@@H](O[C@H]5O[C@H](COS([O-])(=O)=O)[C@@H](O)[C@H](O)[C@H]5NS([O-])(=O)=O)[C@H](O)[C@H]4O)C([O-])=O)[C@H](OS([O-])(=O)=O)[C@H]3NS([O-])(=O)=O)[C@H](O)[C@H]2OS([O-])(=O)=O)C([O-])=O)[C@H](O)[C@H]1NS([O-])(=O)=O',
                                      'name': 'fondaparinux(10-)',
                                      'reason': 'Multiple (2) carboxylate '
                                                'groups found'},
                                  {   'smiles': '[C@H]1([C@H]([C@H](CC1=O)SC[C@H](NC(CC[C@H]([NH3+])C(=O)[O-])=O)C(=O)NCC(=O)[O-])/C=C/[C@H](CCCCC)O)C/C=C\\CCCC(=O)[O-]',
                                      'name': '(S)-PGA2-S-glutathione '
                                              'conjugate(2-)',
                                      'reason': 'Multiple (3) carboxylate '
                                                'groups found'},
                                  {   'smiles': 'C(\\CC=CC=CC=C[C@@H]([C@@H](C/C=C\\CC)O)SC[C@H](NC(CC[C@H]([NH3+])C(=O)[O-])=O)C(=O)NCC(=O)[O-])=C\\C/C=C\\CCC([O-])=O',
                                      'name': '16(S)-glutathionyl-17(R)-hydroxy-(4Z,7Z,10,12,14,19Z)-docosahexaenoate(2-)',
                                      'reason': 'Multiple (3) carboxylate '
                                                'groups found'},
                                  {   'smiles': 'O[C@@H]1[C@H](C(O[C@@H]([C@@H]1O)COS([O-])(=O)=O)O)NC(C)=O',
                                      'name': 'N-acetyl-D-6-sulfogalactosamine(4-)',
                                      'reason': 'No carboxylate group found'},
                                  {   'smiles': 'C1[C@@H]([C@H]([C@H]([C@@H]1O)C/C=C\\CCCCC)/C=C/[C@@H](CCCC([O-])=O)O[C@@H]2O[C@H](C([O-])=O)[C@H]([C@@H]([C@H]2O)O)O)O',
                                      'name': '5-epi-5-F2t-IsoP-glucuronide(2-)',
                                      'reason': 'Multiple (2) carboxylate '
                                                'groups found'},
                                  {   'smiles': 'CCCCCCCC\\C=C/CCCCCCCCCCCCC(O)C(O)=O',
                                      'name': '2-hydroxynervonic acid',
                                      'reason': 'No carboxylate group found'},
                                  {   'smiles': 'C1[C@H]([NH2+]C(C(C1)(C)C)C(=O)[O-])C(=O)[O-]',
                                      'name': '(6S)-3,3-dimethylpiperidinium-2,6-dicarboxylate',
                                      'reason': 'Multiple (2) carboxylate '
                                                'groups found'}],
    'attempt': 0,
    'reasoning': 'To classify a molecule as a "monocarboxylic acid anion", we '
                 'need to determine if the structure contains a deprotonated '
                 'carboxyl group from a monocarboxylic acid. \n'
                 '\n'
                 'Here’s the step-by-step plan to build the function:\n'
                 '\n'
                 '1. Import necessary modules from RDKit.\n'
                 '2. Define the function `is_monocarboxylic_acid_anion` which '
                 'takes a SMILES string as input.\n'
                 '3. Parse the SMILES string into an RDKit molecule object.\n'
                 '4. Check for the presence of the carboxylate group pattern '
                 '`[O-]C(=O)`, which represents the deprotonated carboxylic '
                 'acid.\n'
                 '5. Ensure that this pattern occurs only once, signifying '
                 'that it comes from a monocarboxylic acid (only one '
                 'carboxylic acid group).\n'
                 '6. Return True with an appropriate reason if the conditions '
                 'are met, otherwise return False with the corresponding '
                 'reason.\n'
                 '\n'
                 "Let's implement this in code.",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 1884,
    'num_false_positives': 100,
    'num_true_negatives': 6936,
    'num_false_negatives': 49,
    'num_negatives': None,
    'precision': 0.9495967741935484,
    'recall': 0.9746508018623901,
    'f1': 0.9619606841970897,
    'accuracy': 0.9833872226558145,
    'negative_predictive_value': 0.9929849677881174}