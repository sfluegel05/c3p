"""
Classifies: CHEBI:32952 amine
"""
from rdkit import Chem

def is_amine(smiles: str):
    """
    Determines if a molecule is an amine based on its SMILES string.
    An amine contains a nitrogen atom with at least one alkyl or aryl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define amine SMARTS patterns
    primary_amine_pattern = Chem.MolFromSmarts("[NH2]")
    secondary_amine_pattern = Chem.MolFromSmarts("[NH]([C,R])")
    tertiary_amine_pattern = Chem.MolFromSmarts("[N]([C,R])([C,R])")
    
    # Check for amine groups
    if mol.HasSubstructMatch(primary_amine_pattern):
        return True, "Contains primary amine group (NH2)"
    if mol.HasSubstructMatch(secondary_amine_pattern):
        return True, "Contains secondary amine group (NH)"
    if mol.HasSubstructMatch(tertiary_amine_pattern):
        return True, "Contains tertiary amine group (N)"
    
    return False, "No amine functional group found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:32952',
                          'name': 'amine',
                          'definition': 'A compound formally derived from '
                                        'ammonia by replacing one, two or '
                                        'three hydrogen atoms by hydrocarbyl '
                                        'groups.',
                          'parents': ['CHEBI:50047'],
                          'xrefs': ['KEGG:C00706'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'O1[C@]2([C@@H](O)[C@H](O)[C@@H](O)[C@@]1(OCC[C@H](CCC=C(C(OC[C@]3(O)[C@@H](O)[C@@](OC3)(OC2)[H])=O)C)C)[H])[H]',
                                     'name': 'Urceolide',
                                     'reason': 'No amine functional group '
                                               'found'},
                                 {   'smiles': 'CC1=CC(=C(C=C1)C)C(=O)CSC2=NN=C(S2)C',
                                     'name': '1-(2,5-dimethylphenyl)-2-[(5-methyl-1,3,4-thiadiazol-2-yl)thio]ethanone',
                                     'reason': 'No amine functional group '
                                               'found'},
                                 {   'smiles': 'O=C1C=C([C@@]2(C(C[C@@H](C2)O)(C)C)C)CC[C@@]1(O)C',
                                     'name': 'Enokipodin H',
                                     'reason': 'No amine functional group '
                                               'found'},
                                 {   'smiles': 'CCCC[C@](Cn1cncn1)(C#N)c1ccc(Cl)cc1',
                                     'name': '(S)-myclobutanil',
                                     'reason': 'No amine functional group '
                                               'found'},
                                 {   'smiles': 'OC(=O)CCCCCCC#CCCCC',
                                     'name': '8-tridecynoic acid',
                                     'reason': 'No amine functional group '
                                               'found'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@H](O[C@@H]3[C@@H](CO)O[C@H](O[C@@H]4[C@@H](CO)O[C@H](O)[C@H](O)[C@H]4O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O',
                                     'name': 'alpha-maltotetraose',
                                     'reason': 'No amine functional group '
                                               'found'},
                                 {   'smiles': 'S(C(C)C(O)=O)CC1=CC=CC=C1',
                                     'name': '2-(Benzylthio)propanoic acid',
                                     'reason': 'No amine functional group '
                                               'found'},
                                 {   'smiles': 'OC1=CC=C2C3(OC(C4=C3C=C(C=C4)C(O)=O)=O)C=5C(OC2=C1)=CC(=CC5)O.OC1=CC=C2C3(OC(C4=C3C=CC(=C4)C(O)=O)=O)C=5C(OC2=C1)=CC(=CC5)O',
                                     'name': '5(6)-carboxyfluorescein',
                                     'reason': 'No amine functional group '
                                               'found'},
                                 {   'smiles': '[H][C@@]1(CC[C@@]2([H])[C@]3([H])CC(=O)[C@@]4([H])C[C@H](O)[C@H](O)C[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)[C@@H](O)[C@H](O)[C@@H](C)C(C)C',
                                     'name': 'castasterone',
                                     'reason': 'No amine functional group '
                                               'found'},
                                 {   'smiles': '[H][C@]12CC[C@@](C)(CC1=CC[C@@]1([H])[C@](C)(CO)CCC[C@]21C)C=C',
                                     'name': 'isopimara-7,15-dienol',
                                     'reason': 'No amine functional group '
                                               'found'}],
    'sample_false_negatives': [   {   'smiles': 'CC(=O)OC[C@H]([C@@H](OC(C)=O)c1ccccc1)[N+]([O-])=O',
                                      'name': 'Fenitropan',
                                      'reason': 'No amine functional group '
                                                'found'},
                                  {   'smiles': 'C(=C\\SC)(\\C1=C(OCCCCC)C=CC=C1)/N2C=CN=C2',
                                      'name': 'neticonazole',
                                      'reason': 'No amine functional group '
                                                'found'},
                                  {   'smiles': 'C[C@@H](c1nc(cs1)-c1ccc(cc1)C#N)[C@@](Cn1cncn1)(OCOP(O)(O)=O)c1ccc(F)cc1F',
                                      'name': 'BMS 379224',
                                      'reason': 'No amine functional group '
                                                'found'},
                                  {   'smiles': 'C/C(=C\\CC\\C(\\C)=C\\C[C@]12C(C=C(C([C@H]1O2)OC(C)=O)C)=O)/CCCC(C)(C)O',
                                      'name': 'yanuthone M',
                                      'reason': 'No amine functional group '
                                                'found'},
                                  {   'smiles': 'c1cc[nH]c1',
                                      'name': '1H-pyrrole',
                                      'reason': 'No amine functional group '
                                                'found'},
                                  {   'smiles': 'COc1cc(CCC(=O)n2cccc2)cc(OC)c1OC',
                                      'name': '3-(3,4,5-timethoxyphenyl)propanoylpyrrole',
                                      'reason': 'No amine functional group '
                                                'found'},
                                  {   'smiles': 'OC(Cn1cncn1)(Cn1cncn1)c1ccc(Cl)cc1Cl',
                                      'name': 'UK-47265',
                                      'reason': 'No amine functional group '
                                                'found'},
                                  {   'smiles': 'ClC1=CC=C(C[C@H](N2N=CN=C2)C(O)C(C)(C)C)C=C1',
                                      'name': '(2S)-1-(4-Chlorophenyl)-4,4-dimethyl-2-(1,2,4-triazol-1-yl)-3-pentanol',
                                      'reason': 'No amine functional group '
                                                'found'},
                                  {   'smiles': 'C1=CC=CC=C1CC(/N=C/C(Cl)(Cl)Cl)C',
                                      'name': 'amfecloral',
                                      'reason': 'No amine functional group '
                                                'found'}],
    'attempt': 0,
    'reasoning': 'To classify whether a given chemical entity is an amine '
                 'based on its SMILES string, we need to identify the presence '
                 'of an amine functional group in the molecule. An amine is '
                 'characterized by the presence of a nitrogen atom bonded to '
                 'one or more alkyl or aryl groups. Amines can be primary, '
                 'secondary, or tertiary depending on whether the nitrogen '
                 'atom is bonded to one, two, or three carbon atoms, '
                 'respectively.\n'
                 '\n'
                 'Steps to classify a SMILES string as an amine:\n'
                 '1. Parse the SMILES string to obtain the molecular '
                 'structure.\n'
                 '2. Identify the presence of nitrogen atoms in the '
                 'structure.\n'
                 '3. Check if the nitrogen atoms are part of an amine group:\n'
                 '   - Primary amine: NH2 group\n'
                 '   - Secondary amine: NH group bonded to a carbon atom\n'
                 '   - Tertiary amine: N bonded to three carbon atoms\n'
                 '4. Return True if an amine group is found; otherwise, return '
                 'False.\n'
                 '\n'
                 'The function `is_amine` will utilize RDKit to parse the '
                 'SMILES and identify the amine functionalities.',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 672,
    'num_false_positives': 100,
    'num_true_negatives': 70,
    'num_false_negatives': 9,
    'num_negatives': None,
    'precision': 0.8704663212435233,
    'recall': 0.986784140969163,
    'f1': 0.9249827942188575,
    'accuracy': 0.8719153936545241,
    'negative_predictive_value': 0.8860759493670886}