"""
Classifies: CHEBI:50477 butyrate ester
"""
from rdkit import Chem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule is a butyrate ester based on its SMILES string.
    A butyrate ester contains a butyric acid moiety ('OC(=O)CCC') as part of an ester group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a butyrate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for butyrate ester group
    butyrate_pattern = Chem.MolFromSmarts("OC(=O)CCC")
    
    # Check if molecule has the butyrate ester substructure
    if mol.HasSubstructMatch(butyrate_pattern):
        return True, "Contains butyrate ester substructure"
    else:
        return False, "Does not contain butyrate ester substructure"

# Example usage:
# result, reason = is_butyrate_ester("O=C(OCCCCC)CCC")  # Example, replace with actual SMILES
# print(result, reason)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:50477',
                          'name': 'butyrate ester',
                          'definition': 'Any carboxylic ester where the '
                                        'carboxylic acid component is butyric '
                                        'acid.',
                          'parents': ['CHEBI:33308', 'CHEBI:35748'],
                          'xrefs': ['PMID:23383323'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'O([C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)CO)[C@H]3[C@@H](O)[C@H](OC(O)[C@@H]3O)CO',
                                     'name': 'beta-D-Glcp-(1->2)-beta-D-Glcp-(1->3)-D-Galp',
                                     'reason': 'Does not contain butyrate '
                                               'ester substructure'},
                                 {   'smiles': 'CCS(=O)(=O)NCC[C@@H]1CC[C@@H]([C@H](O1)CO)NC(=O)NC2=CC(=CC(=C2)Cl)Cl',
                                     'name': '1-(3,5-dichlorophenyl)-3-[(2S,3S,6S)-6-[2-(ethylsulfonylamino)ethyl]-2-(hydroxymethyl)-3-oxanyl]urea',
                                     'reason': 'Does not contain butyrate '
                                               'ester substructure'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@H]2[C@@H](CO)O[C@@H](O[C@@H]3[C@@H](CO)O[C@@H](O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O',
                                     'name': 'alpha-D-Galp-(1->4)-beta-D-Galp-(1->4)-beta-D-Glcp',
                                     'reason': 'Does not contain butyrate '
                                               'ester substructure'},
                                 {   'smiles': '[Li+].[Br-]',
                                     'name': 'lithium bromide',
                                     'reason': 'Does not contain butyrate '
                                               'ester substructure'},
                                 {   'smiles': 'C=1(OC)C2=C(C=C3C[C@H]([C@](CC=4C=C(OC)C(OC)=C(C4C13)OC)(C)O)C)OCO2',
                                     'name': 'Besigomsin',
                                     'reason': 'Does not contain butyrate '
                                               'ester substructure'},
                                 {   'smiles': 'C1=CC=C(C(=C1)C=CC(=O)C2=CC=CN2)Cl',
                                     'name': '3-(2-chlorophenyl)-1-(1H-pyrrol-2-yl)-2-propen-1-one',
                                     'reason': 'Does not contain butyrate '
                                               'ester substructure'},
                                 {   'smiles': '[H][C@@]1(COC(C)(C)[C@@]1([H])OC(=O)\\C=C\\c1ccccc1)c1c(O)ccc2C(=O)C[C@H](Oc12)c1ccccc1',
                                     'name': '(+)-tephrorin B',
                                     'reason': 'Does not contain butyrate '
                                               'ester substructure'},
                                 {   'smiles': 'C[C@@H]1CN(C(=O)C2=C(C3=CC=CC=C3CO[C@@H]1CN(C)CC4=CC=CC(=C4)C(=O)O)C5=CC=CC=C5N2C)[C@@H](C)CO',
                                     'name': 'LSM-9341',
                                     'reason': 'Does not contain butyrate '
                                               'ester substructure'},
                                 {   'smiles': 'FC(F)(F)c1nc2c(Cl)c(Cl)c(Cl)c(Cl)c2[nH]1',
                                     'name': '4,5,6,7-Tetrachloro-2-trifluoromethylbenzimidazole',
                                     'reason': 'Does not contain butyrate '
                                               'ester substructure'},
                                 {   'smiles': 'Cl.OC[C@H]1O[C@@H]2[C@@H](Oc3nc(=N)ccn23)[C@@H]1O',
                                     'name': 'ancitabine hydrochloride',
                                     'reason': 'Does not contain butyrate '
                                               'ester substructure'}],
    'sample_false_negatives': [],
    'attempt': 0,
    'reasoning': 'In order to classify a molecule as a butyrate ester based on '
                 'its SMILES string, we need to verify that the SMILES string '
                 'denotes a molecule containing a butyric acid ester, where '
                 'butyric acid is a 4-carbon saturated carboxylic acid '
                 '(butanoic acid). Thus, the molecule should contain the '
                 'substructure –OC(=O)CCC, which represents the ester linkage '
                 'stemming from butyric acid. \n'
                 '\n'
                 "Let's follow these steps for our solution:\n"
                 '\n'
                 '1. Import necessary modules from RDKit to handle molecular '
                 'structures.\n'
                 '2. Define the function `is_butyrate_ester` which takes a '
                 'SMILES string as an argument.\n'
                 '3. Convert the SMILES string to a molecular object.\n'
                 '4. Define the SMARTS pattern that represents the butyrate '
                 'ester substructure.\n'
                 '5. Check if the molecule contains the butyrate ester '
                 'pattern.\n'
                 '6. Return a boolean indicating whether it is a butyrate '
                 'ester along with a reason for this classification.\n'
                 '\n'
                 'Here is the complete code doing the above:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 47,
    'num_false_positives': 100,
    'num_true_negatives': 270,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.3197278911564626,
    'recall': 1.0,
    'f1': 0.4845360824742268,
    'accuracy': 0.7601918465227818,
    'negative_predictive_value': 1.0}