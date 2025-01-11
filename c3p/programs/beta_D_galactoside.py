"""
Classifies: CHEBI:28034 beta-D-galactoside
"""
"""
Classifies: CHEBI:18133 beta-D-galactoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_D_galactoside(smiles: str):
    """
    Determines if a molecule is a beta-D-galactoside based on its SMILES string.
    A beta-D-galactoside is a galactose ring in the beta-configuration, where the hydroxyl group at the anomeric carbon (C1) is in the axial position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-galactoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the beta-D-galactoside pattern
    # The pattern matches a galactose ring with the beta-configuration at the anomeric center
    beta_D_galactoside_pattern = Chem.MolFromSmarts("[C@@H]1([C@@H]([C@H]([C@H]([C@H](O1)CO)O)O)O)")
    
    # Check if the molecule contains the beta-D-galactoside pattern
    if mol.HasSubstructMatch(beta_D_galactoside_pattern):
        return True, "Contains a galactose ring in the beta-configuration"
    else:
        return False, "No galactose ring in the beta-configuration found"

# Example usage:
# smiles = "CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O"  # Methyl beta-D-galactoside
# print(is_beta_D_galactoside(smiles))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:28034',
                          'name': 'beta-D-galactoside',
                          'definition': 'Any D-galactoside having '
                                        'beta-configuration at its anomeric '
                                        'centre.',
                          'parents': ['CHEBI:20961'],
                          'xrefs': ['KEGG:C00602'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'OC[C@H]1O[C@H](C[C@@H]1O)N1C=NC2=C1N=CNC[C@H]2O',
                                     'name': 'pentostatin',
                                     'reason': 'No galactose ring in the '
                                               'beta-configuration found'},
                                 {   'smiles': 'O=C1O[C@@H]([C@@H](OC)C=CC=C(C[C@@H](C)[C@@H]([C@@H]([C@@H]([C@@H](C=C(C=C1OC)C)C)O)C)O)C)[C@H]([C@@H](O)[C@@H]([C@@]2(O[C@H](/C=C/C)[C@@H](C)[C@@H](C2)OC3OC(C(O)C(C3)O)C)O)C)C',
                                     'name': 'Concanamycin D',
                                     'reason': 'No galactose ring in the '
                                               'beta-configuration found'},
                                 {   'smiles': 'O=C(O)[C@]1([C@H]2[C@@](OC=3C=C4C5=C(C=CC=C5)NC4=CC3CC2)(CC[C@@H]1O)C)C',
                                     'name': 'Oxiamycin',
                                     'reason': 'No galactose ring in the '
                                               'beta-configuration found'},
                                 {   'smiles': 'C1CCC(C1)CC#CC2=CC=C(C=C2)[C@@H]3[C@@H]4CN(CC(=O)N4[C@@H]3CO)C(=O)C5=CC=C(C=C5)F',
                                     'name': '(6R,7R,8S)-7-[4-(3-cyclopentylprop-1-ynyl)phenyl]-4-[(4-fluorophenyl)-oxomethyl]-8-(hydroxymethyl)-1,4-diazabicyclo[4.2.0]octan-2-one',
                                     'reason': 'No galactose ring in the '
                                               'beta-configuration found'},
                                 {   'smiles': 'S(=O)(C(SSCCC)CC)CCC',
                                     'name': 'Propyl 1-(propylsulfinyl)propyl '
                                             'disulfide',
                                     'reason': 'No galactose ring in the '
                                               'beta-configuration found'},
                                 {   'smiles': 'ClC=1C(=O)[C@@H]([C@@](O)(C/C=C\\CCCCC)C1)C[C@H](OC(=O)C)[C@@H](OC(=O)C)CCCC(OC)=O',
                                     'name': 'punaglandin 6',
                                     'reason': 'No galactose ring in the '
                                               'beta-configuration found'},
                                 {   'smiles': 'CCCCCCCCCCCCCCCCCCCCCCCC[C@H](O)C(=O)N[C@@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)[C@H](O)[C@@H](O)CCCCCCCCCCCCCC',
                                     'name': 'Ins-1-P-Cer(t18:0/2-OH-26:0)',
                                     'reason': 'No galactose ring in the '
                                               'beta-configuration found'},
                                 {   'smiles': 'CCC(C)(C)C(=O)OC1CC(C=C2C1[C@H]([C@H](C=C2)C)CC[C@@H]3CC(CC(=O)O3)O)C',
                                     'name': '2,2-dimethylbutanoic acid '
                                             '[(7S,8S)-8-[2-[(2R)-4-hydroxy-6-oxo-2-oxanyl]ethyl]-3,7-dimethyl-1,2,3,7,8,8a-hexahydronaphthalen-1-yl] '
                                             'ester',
                                     'reason': 'No galactose ring in the '
                                               'beta-configuration found'},
                                 {   'smiles': '[C@@H]1(C2=CC=CC=C2)[C@@H](N)C1',
                                     'name': '(1S,2R)-tranylcypromine',
                                     'reason': 'No galactose ring in the '
                                               'beta-configuration found'},
                                 {   'smiles': 'OC(=O)/C=C/C#CC#CC#N',
                                     'name': '(2E)-7-Cyanohept-2-en-4,6-diynoic '
                                             'acid',
                                     'reason': 'No galactose ring in the '
                                               'beta-configuration found'}],
    'sample_false_negatives': [],
    'attempt': 0,
    'reasoning': 'To classify a molecule as a beta-D-galactoside, we need to '
                 'identify the presence of a galactose ring in the '
                 'beta-configuration. The beta-configuration means that the '
                 'hydroxyl group at the anomeric carbon (C1) is in the axial '
                 'position, which is typically represented in SMILES as '
                 '`[C@@H]` or `[C@H]` depending on the stereochemistry.\n'
                 '\n'
                 'Here’s the Python program that classifies a molecule as a '
                 'beta-D-galactoside:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 46,
    'num_false_positives': 100,
    'num_true_negatives': 751,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.3150684931506849,
    'recall': 1.0,
    'f1': 0.4791666666666667,
    'accuracy': 0.8885172798216276,
    'negative_predictive_value': 1.0}