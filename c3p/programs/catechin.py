"""
Classifies: CHEBI:23053 catechin
"""
"""
Classifies: catechin
"""
from rdkit import Chem

def is_catechin(smiles: str):
    """
    Determines if a molecule is a catechin based on its SMILES string.
    A catechin is a hydroxyflavan with a flavan-3-ol skeleton and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catechin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES to molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for flavan-3-ol skeleton without stereochemistry
    # This pattern represents the core flavan-3-ol structure
    flavan_3_ol_smarts = 'c1ccc2c(c1)CC(O)CO2'
    flavan_3_ol_pattern = Chem.MolFromSmarts(flavan_3_ol_smarts)
    if flavan_3_ol_pattern is None:
        return False, "Invalid SMARTS pattern for flavan-3-ol skeleton"
    
    # Check if molecule has the flavan-3-ol substructure
    if not mol.HasSubstructMatch(flavan_3_ol_pattern):
        return False, "Flavan-3-ol skeleton not found"
    
    return True, "Molecule contains flavan-3-ol skeleton characteristic of catechins"

# Example usage:
# smiles_example = 'O[C@H]1Cc2c(O)cc(O)cc2O[C@@H]1c1ccc(O)c(O)c1'  # (+)-catechin
# result = is_catechin(smiles_example)
# print(result)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:23053',
                          'name': 'catechin',
                          'definition': 'Members of the class of hydroxyflavan '
                                        'that have a flavan-3-ol skeleton and '
                                        'its substituted derivatives.',
                          'parents': ['CHEBI:72010'],
                          'xrefs': ['KEGG:C17590', 'LINCS:LSM-1682'],
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
               'O1C(C2C(OC3OC(C(O)C(O)C3O)CO)C1(OC=4C2=C(O)C=C(O)C4)C=5C=C(OC)C(O)=C(O)C5)(C)C(=O)C '
               "NAME: 5'-Hydroxycastavinol REASON: MISSED Flavan-3-ol skeleton "
               'not found\n'
               ' * SMILES: '
               'S([C@@H]1[C@@H](OC(=O)C2=CC(O)=C(O)C(O)=C2)[C@H](OC=3C1=C(O)C=C(O)C3)C4=CC(O)=C(O)C=C4)CCN '
               'NAME: 4beta-(2-Aminoethylthio)epicatechin 3-gallate REASON: '
               'MISSED Flavan-3-ol skeleton not found\n'
               ' * SMILES: O[C@H]1Cc2c(O)cc(O)cc2O[C@@H]1c1ccc(O)c(O)c1 NAME: '
               '(+)-catechin REASON: MISSED Flavan-3-ol skeleton not found\n'
               ' * SMILES: '
               'O[C@H]1Cc2c(O)cc(O)c([C@H]3[C@H](O)[C@H](Oc4cc(O)ccc34)c3cc(O)c(O)c(O)c3)c2O[C@@H]1c1cc(O)c(O)c(O)c1 '
               'NAME: robinetinidol-(4alpha,8)-gallocatechin REASON: MISSED '
               'Flavan-3-ol skeleton not found\n'
               ' * SMILES: '
               'COc1ccc(cc1O)[C@H]1Oc2cc(O)cc(O)c2C[C@H]1OC(=O)c1cc(O)c(O)c(O)c1 '
               "NAME: 4'-O-methylepicatechin-3-O-gallate REASON: MISSED "
               'Flavan-3-ol skeleton not found\n'
               ' * SMILES: '
               'O[C@@H](Cc1c(O)cc(O)cc1O)[C@H](c1ccc(O)c(O)c1)c1c(O)cc2O[C@@H]([C@@H](O)Cc2c1O)c1ccc(O)c(O)c1 '
               'NAME: Gambiriin A3 REASON: MISSED Flavan-3-ol skeleton not '
               'found\n'
               ' * SMILES: '
               'O1C(C(O)CC=2C1=CC(O)=CC2O)C=3C=4C(C(O)=C(O)C3)=C(O)C(=O)C=C(C4)C(O)=O '
               'NAME: Theaflavic acid REASON: MISSED Flavan-3-ol skeleton not '
               'found\n'
               ' * SMILES: '
               'S(OC1=CC(C2OC=3C(C(=O)C2OC(=O)C4=CC(O)=C(O)C(O)=C4)=C(O)C=C(O)C3)=CC(O)=C1O)(O)(=O)=O '
               'NAME: Myricatin REASON: MISSED Flavan-3-ol skeleton not found\n'
               ' * SMILES: '
               'O[C@H]1Cc2c(O)cc3O[C@@H](Cc4c(O)cc(O)cc4O)[C@@H](c3c2O[C@@H]1c1ccc(O)c(O)c1)c1ccc(O)c(O)c1 '
               'NAME: Gambiriin B1 REASON: MISSED Flavan-3-ol skeleton not '
               'found\n'
               ' * SMILES: '
               'S(OC=1C=C2O[C@@H]([C@H](O)CC2=C(O)C1)C3=CC(OC)=C(O)C=C3)(O)(=O)=O '
               "NAME: 3'-methylepicatechin-7-sulfate REASON: MISSED "
               'Flavan-3-ol skeleton not found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C[C@@H]1CN(S(=O)(=O)C2=C(C=C(C=C2)C3=CC=C(C=C3)C(=O)N(C)C)O[C@@H]1CN(C)C(=O)NC4=CC=CC=C4F)[C@H](C)CO',
                                     'name': '4-[(4R,5S)-5-[[[(2-fluoroanilino)-oxomethyl]-methylamino]methyl]-2-[(2R)-1-hydroxypropan-2-yl]-4-methyl-1,1-dioxo-4,5-dihydro-3H-6,1$l^{6},2-benzoxathiazocin-8-yl]-N,N-dimethylbenzamide',
                                     'reason': 'Flavan-3-ol skeleton not '
                                               'found'},
                                 {   'smiles': 'C[C@H]1[C@H](CN(C(=O)C2=C(C=CC(=C2)C#N)OC[C@H]3[C@H](CC[C@@H](O3)CCN(C1=O)C)OC)C)OC',
                                     'name': 'LSM-10564',
                                     'reason': 'Flavan-3-ol skeleton not '
                                               'found'},
                                 {   'smiles': 'CCC(C(=O)NC1=NN=C(S1)SCC)(C(F)(F)F)C(F)(F)F',
                                     'name': 'N-[5-(ethylthio)-1,3,4-thiadiazol-2-yl]-2,2-bis(trifluoromethyl)butanamide',
                                     'reason': 'Flavan-3-ol skeleton not '
                                               'found'},
                                 {   'smiles': 'O=C(N[C@@H](CCC(=O)N)C(O)=O)[C@@H](NC(=O)[C@@H](N)CC1=CC=C(O)C=C1)CC(=O)N',
                                     'name': 'Tyr-Asn-Gln',
                                     'reason': 'Flavan-3-ol skeleton not '
                                               'found'},
                                 {   'smiles': 'OC[C@H]1OC(O)[C@@H](O)[C@@H]1O',
                                     'name': 'D-arabinofuranose',
                                     'reason': 'Flavan-3-ol skeleton not '
                                               'found'},
                                 {   'smiles': 'O[C@@H]([C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\\3/C[C@@H](O)CCC3=C)[H])C)[H])C)C',
                                     'name': '(5Z,7E)-(3S,22R)- '
                                             '24-nor-9,10-seco-5,7,10(19)-cholatriene-3,22-diol',
                                     'reason': 'Flavan-3-ol skeleton not '
                                               'found'},
                                 {   'smiles': 'O1C(C(CC1=O)C)CCCCCC',
                                     'name': '5-Hexyldihydro-4-methyl-2(3H)-furanone',
                                     'reason': 'Flavan-3-ol skeleton not '
                                               'found'},
                                 {   'smiles': 'O=C1O[C@]23C(=O)N[C@H]([C@@H]2[C@H](C)C(=C[C@@H]3C=CC[C@@H]([C@H]([C@H](C=C1)C)O)C)C)CC4=CC=CC=C4',
                                     'name': '7-deoxy-cytochalasin Z7',
                                     'reason': 'Flavan-3-ol skeleton not '
                                               'found'},
                                 {   'smiles': '[C@]123[C@H](O1)[C@](C(C=C2[C@]4([C@](CC3)(C([C@@H](C(C4)=O)O)(C)C)[H])C)=O)(C=C)C',
                                     'name': 'yucalexin P8',
                                     'reason': 'Flavan-3-ol skeleton not '
                                               'found'},
                                 {   'smiles': '[H][C@]12CC[C@@]3([H])[C@]4([H])CC[C@]([H])([C@H](C)CCC(O)=O)[C@@]4(C)[C@@H](O)C[C@]3([H])[C@@]1(C)CC[C@@H](O)C2',
                                     'name': 'deoxycholic acid',
                                     'reason': 'Flavan-3-ol skeleton not '
                                               'found'}],
    'sample_false_negatives': [],
    'attempt': 1,
    'reasoning': 'The previous program failed because the SMARTS pattern used '
                 'to identify the flavan-3-ol skeleton was incorrect or too '
                 'restrictive. As a result, it did not match the structures of '
                 'known catechins provided. The SMARTS pattern may have been '
                 'too specific or improperly formatted, leading to no '
                 'matches.\n'
                 '\n'
                 'To improve the program, we need to define a more accurate '
                 'and general SMARTS pattern that correctly represents the '
                 'flavan-3-ol skeleton found in catechins. By creating a new '
                 'SMARTS pattern that captures the core structure of catechins '
                 'without enforcing stereochemistry constraints, we can ensure '
                 'it matches catechin molecules and their substituted '
                 'derivatives. We can test this new pattern on the example '
                 'catechin SMILES strings to verify that it works as expected.',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 106,
    'num_false_positives': 100,
    'num_true_negatives': 27714,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.5145631067961165,
    'recall': 1.0,
    'f1': 0.6794871794871795,
    'accuracy': 0.9964183381088825,
    'negative_predictive_value': 1.0}