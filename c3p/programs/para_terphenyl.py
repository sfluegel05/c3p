"""
Classifies: CHEBI:75874 para-terphenyl
"""
from rdkit import Chem

def is_para_terphenyl(smiles: str):
    """
    Determines if a molecule is a para-terphenyl based on its SMILES string.
    A para-terphenyl consists of a 1,4-diphenylbenzene skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a para-terphenyl, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for para-terphenyl: c1ccc(cc1)-c1ccc(cc1)-c1ccccc1
    para_terphenyl_pattern = Chem.MolFromSmarts("c1ccc(cc1)-c1ccc(cc1)-c1ccccc1")
    
    # Check if the molecule matches the para-terphenyl pattern
    if mol.HasSubstructMatch(para_terphenyl_pattern):
        return True, "Contains 1,4-diphenylbenzene skeleton"
    else:
        return False, "Does not contain 1,4-diphenylbenzene skeleton"

# Example usage
example_smiles = "c1ccc(cc1)-c1ccc(cc1)-c1ccccc1"  # 1,4-diphenylbenzene itself
result, reason = is_para_terphenyl(example_smiles)
print(result, reason)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:75874',
                          'name': 'para-terphenyl',
                          'definition': 'A ring assembly based on a '
                                        '1,4-diphenylbenzene skeleton and its '
                                        'substituted derivatives thereof.',
                          'parents': ['CHEBI:33836', 'CHEBI:36820'],
                          'xrefs': ['Wikipedia:Terphenyl'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'O([C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)CO)[C@H]3[C@@H](O)[C@H](OC(O)[C@@H]3O)CO',
                                     'name': 'beta-D-Glcp-(1->2)-beta-D-Glcp-(1->3)-D-Galp',
                                     'reason': 'Does not contain '
                                               '1,4-diphenylbenzene skeleton'},
                                 {   'smiles': 'CCS(=O)(=O)NCC[C@@H]1CC[C@@H]([C@H](O1)CO)NC(=O)NC2=CC(=CC(=C2)Cl)Cl',
                                     'name': '1-(3,5-dichlorophenyl)-3-[(2S,3S,6S)-6-[2-(ethylsulfonylamino)ethyl]-2-(hydroxymethyl)-3-oxanyl]urea',
                                     'reason': 'Does not contain '
                                               '1,4-diphenylbenzene skeleton'},
                                 {   'smiles': 'C(C(O)=O)C/C=C\\C/C=C\\C\\C=C/C=C/C=C/[C@H]1[C@H](C/C=C\\CC)O1',
                                     'name': '(16S,17S)-epoxy-(4Z,7Z,10Z,12E,14E,19Z)-docosahexaenoic '
                                             'acid',
                                     'reason': 'Does not contain '
                                               '1,4-diphenylbenzene skeleton'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@H]2[C@@H](CO)O[C@@H](O[C@@H]3[C@@H](CO)O[C@@H](O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O',
                                     'name': 'alpha-D-Galp-(1->4)-beta-D-Galp-(1->4)-beta-D-Glcp',
                                     'reason': 'Does not contain '
                                               '1,4-diphenylbenzene skeleton'},
                                 {   'smiles': '[Li+].[Br-]',
                                     'name': 'lithium bromide',
                                     'reason': 'Does not contain '
                                               '1,4-diphenylbenzene skeleton'},
                                 {   'smiles': 'C=1(OC)C2=C(C=C3C[C@H]([C@](CC=4C=C(OC)C(OC)=C(C4C13)OC)(C)O)C)OCO2',
                                     'name': 'Besigomsin',
                                     'reason': 'Does not contain '
                                               '1,4-diphenylbenzene skeleton'},
                                 {   'smiles': 'C1=CC=C(C(=C1)C=CC(=O)C2=CC=CN2)Cl',
                                     'name': '3-(2-chlorophenyl)-1-(1H-pyrrol-2-yl)-2-propen-1-one',
                                     'reason': 'Does not contain '
                                               '1,4-diphenylbenzene skeleton'},
                                 {   'smiles': 'O=C1OC(O)C(=C1C(O)C(C)C)C',
                                     'name': '5-Hydroxy-3-(1-hydroxy-2-methylpropyl)-4-methyl-2(5H)-furanone',
                                     'reason': 'Does not contain '
                                               '1,4-diphenylbenzene skeleton'},
                                 {   'smiles': '[H][C@@]1(COC(C)(C)[C@@]1([H])OC(=O)\\C=C\\c1ccccc1)c1c(O)ccc2C(=O)C[C@H](Oc12)c1ccccc1',
                                     'name': '(+)-tephrorin B',
                                     'reason': 'Does not contain '
                                               '1,4-diphenylbenzene skeleton'},
                                 {   'smiles': 'C[C@@H]1CN(C(=O)C2=C(C3=CC=CC=C3CO[C@@H]1CN(C)CC4=CC=CC(=C4)C(=O)O)C5=CC=CC=C5N2C)[C@@H](C)CO',
                                     'name': 'LSM-9341',
                                     'reason': 'Does not contain '
                                               '1,4-diphenylbenzene skeleton'}],
    'sample_false_negatives': [],
    'attempt': 0,
    'reasoning': 'To create a program that classifies chemical entities of the '
                 'class para-terphenyl, we need a function that can identify '
                 'the 1,4-diphenylbenzene skeleton in a given SMILES string. '
                 'Para-terphenyl and its derivatives have this characteristic '
                 'structure. The classification function `is_para_terphenyl` '
                 'should utilize the RDKit library to identify this core '
                 'structure.\n'
                 '\n'
                 'The steps to create this function are as follows:\n'
                 '\n'
                 '1. **Import the necessary modules from RDKit:** RDKit '
                 'provides cheminformatics functionalities to analyze and '
                 'manipulate molecular structures.\n'
                 '\n'
                 '2. **Define the SMARTS pattern for the para-terphenyl '
                 'structure:** In this case, the pattern for '
                 '1,4-diphenylbenzene can be represented in SMARTS language to '
                 'match three benzene rings connected such that there is a '
                 'central benzene with two phenyl groups at the para (1,4) '
                 'position.\n'
                 '\n'
                 '3. **Create the `is_para_terphenyl` function:** This '
                 'function will:\n'
                 '   - Take a SMILES string as input.\n'
                 '   - Use RDKit to parse the SMILES string and search for the '
                 'para-terphenyl pattern.\n'
                 '   - Return a boolean indicating whether the pattern is '
                 'found, along with a reason for the classification.\n'
                 '\n'
                 '4. **Handle invalid SMILES strings by returning (False, '
                 '"Invalid SMILES string").**\n'
                 '\n'
                 "Let's implement this program:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 74,
    'num_false_positives': 23,
    'num_true_negatives': 142203,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.7628865979381443,
    'recall': 1.0,
    'f1': 0.8654970760233918,
    'accuracy': 0.9998383696416022,
    'negative_predictive_value': 1.0}