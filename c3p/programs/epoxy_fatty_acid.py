"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is an epoxy fatty acid, False otherwise.
        str: Reason for classification.
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the epoxide group pattern
    epoxide_pattern = Chem.MolFromSmarts("[C]1[O][C]1")
    if not mol.HasSubstructMatch(epoxide_pattern):
        return False, "No epoxide group found"

    # Check for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Fatty acids typically have 12-22 carbon atoms
    if c_count < 12 or c_count > 22:
        return False, f"Unexpected carbon count for a fatty acid, found {c_count}"

    # Ensure that the structure is mostly a linear chain with an epoxide
    # Ensure it is not overly cyclic or complex
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings > 1:
        return False, "Unexpected complex ring structure"

    return True, "Identified as an epoxy fatty acid with appropriate functional groups and chain length"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:61498',
                          'name': 'epoxy fatty acid',
                          'definition': 'A heterocyclic fatty acid containing '
                                        'an epoxide ring as part of its '
                                        'structure.',
                          'parents': ['CHEBI:23931', 'CHEBI:48847'],
                          'xrefs': ['PMID:15145985'],
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
               'C(O)(=O)CC/C=C\\C[C@H]1[C@H](/C=C/C=C/C=C\\C=C\\[C@@H](C/C=C\\CC)O)O1 '
               'NAME: 7S,8S-epoxy-17R-hydroxydocosahexaenoic acid REASON: '
               'MISSED Insufficient carbon chain length for fatty acid '
               'structure\n'
               ' * SMILES: '
               'C(C(O)=O)C/C=C\\C/C=C\\C\\C=C/C=C/C=C/[C@H]1[C@H](C/C=C\\CC)O1 '
               'NAME: (16S,17S)-epoxy-(4Z,7Z,10Z,12E,14E,19Z)-docosahexaenoic '
               'acid REASON: MISSED Insufficient carbon chain length for fatty '
               'acid structure\n'
               ' * SMILES: C(CCCO)C/C=C\\C/C=C\\C/C=C\\CC1C(CCCC(O)=O)O1 NAME: '
               '5,6-epoxy-20-hydroxy-(8Z,11Z,14Z)-icosatrienoic acid REASON: '
               'MISSED Insufficient carbon chain length for fatty acid '
               'structure\n'
               ' * SMILES: '
               'O=C(CCC/C=C\\C/C=C\\C/C=C\\[C@@H]([C@H]1[C@H](CCCCC)O1)O)O '
               'NAME: (13S)-hydroxy-(14S,15S)-epoxy-(5Z,8Z,11Z)-icosatrienoic '
               'acid REASON: MISSED Insufficient carbon chain length for fatty '
               'acid structure\n'
               ' * SMILES: '
               'CCCCC[C@@H]1O[C@H]1[C@H](O)\\C=C/C\\C=C/C\\C=C/CCCC(O)=O NAME: '
               '(13R)-hydroxy-(14S,15S)-epoxyicosa-(5Z,8Z,11Z)-trienoic acid '
               'REASON: MISSED Insufficient carbon chain length for fatty acid '
               'structure\n'
               ' * SMILES: '
               'C(CCC(O)=O)/C=C\\C/C=C\\C=C\\[C@H](C[C@H]1[C@@H](CCCCC)O1)OO '
               'NAME: '
               '(12S)-hydroperoxy-(14S,15R)-epoxy-(5Z,8Z,10E)-icosatrienoic '
               'acid REASON: MISSED Insufficient carbon chain length for fatty '
               'acid structure\n'
               ' * SMILES: C1(C(C/C=C\\CCCCCO)O1)CCCCCCCC(=O)O NAME: '
               '9,10-epoxy-18-hydroxy-(12Z)-octadecenoic acid REASON: MISSED '
               'Insufficient carbon chain length for fatty acid structure\n'
               ' * SMILES: O1C(C1C#CC#CCCCCCCCC)CCCC(O)=O NAME: '
               '5,6-Epoxyoctadeca-7,9-diynoic acid REASON: MISSED Insufficient '
               'carbon chain length for fatty acid structure\n'
               ' * SMILES: '
               'CCCCC\\C=C/C[C@@H]1O[C@H]1\\C=C\\[C@H](O)C\\C=C/CCCC(O)=O '
               'NAME: (8R)-hydroxy-(11S,12S)-epoxyicosa-(5Z,9E,14Z)-trienoic '
               'acid REASON: MISSED Insufficient carbon chain length for fatty '
               'acid structure\n'
               ' * SMILES: CC[C@]1(C)O[C@@H]1CC\\C(C)=C\\CC\\C(C)=C\\C(O)=O '
               'NAME: '
               '(2E,6E,10R,11S)-10,11-epoxy-3,7,11-trimethyltrideca-2,6-dienoic '
               'acid REASON: MISSED Insufficient carbon chain length for fatty '
               'acid structure\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'O([C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)CO)[C@H]3[C@@H](O)[C@H](OC(O)[C@@H]3O)CO',
                                     'name': 'beta-D-Glcp-(1->2)-beta-D-Glcp-(1->3)-D-Galp',
                                     'reason': 'No epoxide group found'},
                                 {   'smiles': 'CCS(=O)(=O)NCC[C@@H]1CC[C@@H]([C@H](O1)CO)NC(=O)NC2=CC(=CC(=C2)Cl)Cl',
                                     'name': '1-(3,5-dichlorophenyl)-3-[(2S,3S,6S)-6-[2-(ethylsulfonylamino)ethyl]-2-(hydroxymethyl)-3-oxanyl]urea',
                                     'reason': 'No epoxide group found'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@H]2[C@@H](CO)O[C@@H](O[C@@H]3[C@@H](CO)O[C@@H](O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O',
                                     'name': 'alpha-D-Galp-(1->4)-beta-D-Galp-(1->4)-beta-D-Glcp',
                                     'reason': 'No epoxide group found'},
                                 {   'smiles': '[Li+].[Br-]',
                                     'name': 'lithium bromide',
                                     'reason': 'No epoxide group found'},
                                 {   'smiles': 'C=1(OC)C2=C(C=C3C[C@H]([C@](CC=4C=C(OC)C(OC)=C(C4C13)OC)(C)O)C)OCO2',
                                     'name': 'Besigomsin',
                                     'reason': 'No epoxide group found'},
                                 {   'smiles': 'C1=CC=C(C(=C1)C=CC(=O)C2=CC=CN2)Cl',
                                     'name': '3-(2-chlorophenyl)-1-(1H-pyrrol-2-yl)-2-propen-1-one',
                                     'reason': 'No epoxide group found'},
                                 {   'smiles': 'O=C1OC(O)C(=C1C(O)C(C)C)C',
                                     'name': '5-Hydroxy-3-(1-hydroxy-2-methylpropyl)-4-methyl-2(5H)-furanone',
                                     'reason': 'No epoxide group found'},
                                 {   'smiles': '[H][C@@]1(COC(C)(C)[C@@]1([H])OC(=O)\\C=C\\c1ccccc1)c1c(O)ccc2C(=O)C[C@H](Oc12)c1ccccc1',
                                     'name': '(+)-tephrorin B',
                                     'reason': 'No epoxide group found'},
                                 {   'smiles': 'C[C@@H]1CN(C(=O)C2=C(C3=CC=CC=C3CO[C@@H]1CN(C)CC4=CC=CC(=C4)C(=O)O)C5=CC=CC=C5N2C)[C@@H](C)CO',
                                     'name': 'LSM-9341',
                                     'reason': 'No epoxide group found'},
                                 {   'smiles': 'FC(F)(F)c1nc2c(Cl)c(Cl)c(Cl)c(Cl)c2[nH]1',
                                     'name': '4,5,6,7-Tetrachloro-2-trifluoromethylbenzimidazole',
                                     'reason': 'No epoxide group found'}],
    'sample_false_negatives': [   {   'smiles': 'ClC[C@]1(O)[C@@H]([C@]2(O[C@@H]2CC=C(C)C)C)[C@H](OC)[C@H](OC(=O)CCC(=O)O)CC1',
                                      'name': 'Ligerin',
                                      'reason': 'Unexpected complex ring '
                                                'structure'},
                                  {   'smiles': 'O1[C@H]([C@@H]1[C@H]([C@@H](O)C)C)C[C@@H]2[C@@H](O)[C@@H](O)[C@@H](OC2)C/C(/C)=C/C(O)=O',
                                      'name': 'Monic acid',
                                      'reason': 'Unexpected complex ring '
                                                'structure'}],
    'attempt': 4,
    'reasoning': 'In the previous attempt, it appears that the code failed to '
                 'correctly identify the epoxy fatty acids based on their '
                 'SMILES strings. The key components to focus on are:\n'
                 '\n'
                 '1. **Detection of the Epoxide Ring:** The epoxide group '
                 'structure (a three-membered ether ring) is central to epoxy '
                 'fatty acids. The SMARTS pattern used for epoxide detection '
                 'might have been too strict or incorrect.\n'
                 '\n'
                 '2. **Carbon Chain Verification:** The definition calls for a '
                 '"heterocyclic fatty acid," implying a long hydrocarbon chain '
                 'characteristic of fatty acids. The failure case indicated '
                 'insufficient chain length, which suggests the pattern for '
                 'detecting the carbon chain was either incorrect or too '
                 'restrictive.\n'
                 '\n'
                 '3. **Carboxylic Acid Group:** This group is essential in '
                 'fatty acids and should have a straightforward pattern to '
                 'detect.\n'
                 '\n'
                 "Let's refine the approach:\n"
                 '- Use a correct SMARTS pattern for epoxide rings: '
                 '`"[C]1[O][C]1"`\n'
                 '- Verify the presence of the carboxylic acid group using '
                 '`"[CX3](=O)[OX2H1]"`\n'
                 "- Ensure that the molecule doesn't contain undesired complex "
                 'ring systems beyond the epoxide.\n'
                 '- Count the number of carbons to be within a reasonable '
                 'range for fatty acids (though the previous filtering '
                 'suggests chain length is sometimes misidentified).\n'
                 '\n'
                 "Here's how the updated classifier might look:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 77,
    'num_false_positives': 51,
    'num_true_negatives': 142170,
    'num_false_negatives': 2,
    'num_negatives': None,
    'precision': 0.6015625,
    'recall': 0.9746835443037974,
    'f1': 0.7439613526570048,
    'accuracy': 0.9996275474349965,
    'negative_predictive_value': 0.9999859325324255}