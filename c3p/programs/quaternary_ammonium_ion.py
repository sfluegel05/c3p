"""
Classifies: CHEBI:35267 quaternary ammonium ion
"""
"""
Classifies: quaternary ammonium ion
"""
from rdkit import Chem

def is_quaternary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a quaternary ammonium ion based on its SMILES string.
    A quaternary ammonium ion is a derivative of ammonium, NH4(+), in which all four of the hydrogens bonded to nitrogen have been replaced with univalent (usually organyl) groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quaternary ammonium ion, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define quaternary ammonium ion SMARTS pattern
    # Nitrogen atom with formal charge +1 and degree 4 (four single bonds)
    quat_ammonium_pattern = Chem.MolFromSmarts("[N+;D4]")
    if mol.HasSubstructMatch(quat_ammonium_pattern):
        return True, "Contains quaternary ammonium ion"
    else:
        return False, "No quaternary ammonium ion found"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'quaternary ammonium ion',
        'definition': 'A derivative of ammonium, NH4(+), in which all four of the hydrogens bonded to nitrogen have been replaced with univalent (usually organyl) groups.',
        'parents': []
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': 50,
        'max_negative_in_prompt': 20,
        'max_instances_in_prompt': 100,
        'test_proportion': 0.1
    },
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
}


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35267',
                          'name': 'quaternary ammonium ion',
                          'definition': 'A derivative of ammonium, NH4(+), in '
                                        'which all four of the hydrogens '
                                        'bonded to nitrogen have been replaced '
                                        'with univalent (usually organyl) '
                                        'groups.',
                          'parents': ['CHEBI:25697', 'CHEBI:35274'],
                          'xrefs': ['KEGG:C06703'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': '[H]C(=O)\\C=C\\C=C\\C=C/CC',
                                     'name': '(2E,4E,6Z)-nonatrienal',
                                     'reason': 'No quaternary ammonium ion '
                                               'found'},
                                 {   'smiles': 'C=1(C2=CC=C(Cl)C=C2)NN=C(C1C3=CC=NC=N3)C4CCN(CC4)C(=O)CO',
                                     'name': 'SD-06',
                                     'reason': 'No quaternary ammonium ion '
                                               'found'},
                                 {   'smiles': 'CS(=O)(=O)C1=CC=C(C=C1)C([C@H](CF)NC(=O)C(Cl)Cl)O',
                                     'name': '2,2-dichloro-N-[(2R)-3-fluoro-1-hydroxy-1-(4-methylsulfonylphenyl)propan-2-yl]acetamide',
                                     'reason': 'No quaternary ammonium ion '
                                               'found'},
                                 {   'smiles': 'CC=CC1=CC=C(C=C1)[C@@H]2[C@@H]3CN(CC(=O)N3[C@H]2CO)C(=O)NC4CCCCC4',
                                     'name': 'LSM-38199',
                                     'reason': 'No quaternary ammonium ion '
                                               'found'},
                                 {   'smiles': 'O1C(N2C3=NC=NC(NC4CCCC4)=C3N=C2NC)C(O)C(O)C1CO',
                                     'name': '(2R,5R)-2-[6-(Cyclopentylamino)-8-(methylamino)purin-9-yl]-5-(hydroxymethyl)oxolane-3,4-diol',
                                     'reason': 'No quaternary ammonium ion '
                                               'found'},
                                 {   'smiles': 'SC[C@H](N)C(=O)N[C@@H](CC(O)=O)C(=O)N[C@@H](CO)C(O)=O',
                                     'name': 'Cys-Asp-Ser',
                                     'reason': 'No quaternary ammonium ion '
                                               'found'},
                                 {   'smiles': 'C[C@]12CC[C@H]3[C@H]([C@H]1C[C@H]([C@@H]2O)O)CCC4=C3C=CC(=C4)O',
                                     'name': '(8R,9S,13S,14R,16R,17R)-13-methyl-6,7,8,9,11,12,14,15,16,17-decahydrocyclopenta[a]phenanthrene-3,16,17-triol',
                                     'reason': 'No quaternary ammonium ion '
                                               'found'},
                                 {   'smiles': '[H][C@](COP(OC[C@@](COC(CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)=O)(OC(CCCCCCC/C=C\\CCCCCCCC)=O)[H])(=O)O)(C(O)=O)N',
                                     'name': 'PS(18:3(9Z,12Z,15Z)/18:1(9Z))',
                                     'reason': 'No quaternary ammonium ion '
                                               'found'},
                                 {   'smiles': 'C1CCC(CC1)C(=O)NCC[C@H]2CC[C@@H]([C@@H](O2)CO)NC(=O)NC3=CC=CC=C3F',
                                     'name': 'N-[2-[(2R,5S,6R)-5-[[(2-fluoroanilino)-oxomethyl]amino]-6-(hydroxymethyl)-2-oxanyl]ethyl]cyclohexanecarboxamide',
                                     'reason': 'No quaternary ammonium ion '
                                               'found'},
                                 {   'smiles': 'C1CC(=O)N[C@@H]1C(=O)N2CSC[C@H]2C(=O)O',
                                     'name': '(4R)-3-[oxo-[(2S)-5-oxo-2-pyrrolidinyl]methyl]-4-thiazolidinecarboxylic '
                                             'acid',
                                     'reason': 'No quaternary ammonium ion '
                                               'found'}],
    'sample_false_negatives': [   {   'smiles': 'P(OC[C@H](OC(=O)CCCCCCCCCCCCCCCCCCCC)COC(=O)CCCCCCCCCCCCCCCCCCC)(OC[C@@H](O)CO)(O)=O',
                                      'name': 'PG(20:0/21:0)',
                                      'reason': 'No quaternary ammonium ion '
                                                'found'},
                                  {   'smiles': 'P(OC[C@H](OC(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCCCCCCCC/C=C\\C/C=C\\CCCCC)(O)(O)=O',
                                      'name': 'PA(20:2(11Z,14Z)/22:6(4Z,7Z,10Z,13Z,16Z,19Z))',
                                      'reason': 'No quaternary ammonium ion '
                                                'found'}],
    'attempt': 0,
    'reasoning': 'To determine if a molecule is a quaternary ammonium ion, we '
                 'need to look for a nitrogen atom that is positively charged '
                 'and bonded to four substituents, where all hydrogens have '
                 'been replaced with univalent (usually organic) groups. In '
                 'SMILES notation, quaternary ammonium ions often have the '
                 'nitrogen represented as `[N+]` with four single bonds to '
                 'other atoms/groups.\n'
                 '\n'
                 'We can use RDKit to parse the SMILES string and search for '
                 "nitrogen atoms that meet these criteria. Specifically, we'll "
                 'look for a nitrogen atom with:\n'
                 '\n'
                 '- A formal charge of +1\n'
                 '- A degree of 4 (i.e., bonded to four atoms)\n'
                 '- All bonds are single bonds\n'
                 '\n'
                 'Below is the Python code that implements this logic:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 1329,
    'num_false_positives': 100,
    'num_true_negatives': 19491,
    'num_false_negatives': 2,
    'num_negatives': None,
    'precision': 0.9300209937018894,
    'recall': 0.9984973703981969,
    'f1': 0.9630434782608696,
    'accuracy': 0.9951247490679668,
    'negative_predictive_value': 0.9998973990663315}