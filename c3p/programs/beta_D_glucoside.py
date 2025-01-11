"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
"""
Classifies: CHEBI:59805 beta-D-glucoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_beta_D_glucoside(smiles: str):
    """
    Determines if a molecule is a beta-D-glucoside based on its SMILES string.
    A beta-D-glucoside is a molecule where a beta-D-glucose is attached to another molecule via a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-glucoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the beta-D-glucose pattern
    beta_D_glucose_pattern = Chem.MolFromSmarts("[C@@H]1([C@H]([C@@H]([C@H]([C@@H](O1)CO)O)O)O)O")
    if not mol.HasSubstructMatch(beta_D_glucose_pattern):
        return False, "No beta-D-glucose moiety found"

    # Check if the glucose is attached to another molecule via a glycosidic bond
    glycosidic_bond_pattern = Chem.MolFromSmarts("[C@@H]1([C@H]([C@@H]([C@H]([C@@H](O1)CO)O)O)O)O-[*]")
    if not mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "No glycosidic bond found"

    # Verify the beta configuration by checking the anomeric carbon
    anomeric_carbon_pattern = Chem.MolFromSmarts("[C@@H]1([C@H]([C@@H]([C@H]([C@@H](O1)CO)O)O)O)O")
    matches = mol.GetSubstructMatches(anomeric_carbon_pattern)
    for match in matches:
        anomeric_carbon = mol.GetAtomWithIdx(match[0])
        if anomeric_carbon.GetHybridization() != Chem.HybridizationType.SP3:
            return False, "Anomeric carbon is not in the correct configuration"

    return True, "Contains a beta-D-glucose moiety attached via a glycosidic bond"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:22798',
                          'name': 'beta-D-glucoside',
                          'definition': 'Any D-glucoside in which the anomeric '
                                        'centre has beta-configuration.',
                          'parents': ['CHEBI:35436', 'CHEBI:60980'],
                          'xrefs': ['KEGG:C00963'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'OC[C@H]1O[C@H](C[C@@H]1O)N1C=NC2=C1N=CNC[C@H]2O',
                                     'name': 'pentostatin',
                                     'reason': 'No beta-D-glucose moiety '
                                               'found'},
                                 {   'smiles': 'O=C1O[C@@H]([C@@H](OC)C=CC=C(C[C@@H](C)[C@@H]([C@@H]([C@@H]([C@@H](C=C(C=C1OC)C)C)O)C)O)C)[C@H]([C@@H](O)[C@@H]([C@@]2(O[C@H](/C=C/C)[C@@H](C)[C@@H](C2)OC3OC(C(O)C(C3)O)C)O)C)C',
                                     'name': 'Concanamycin D',
                                     'reason': 'No beta-D-glucose moiety '
                                               'found'},
                                 {   'smiles': 'O=C(O)[C@]1([C@H]2[C@@](OC=3C=C4C5=C(C=CC=C5)NC4=CC3CC2)(CC[C@@H]1O)C)C',
                                     'name': 'Oxiamycin',
                                     'reason': 'No beta-D-glucose moiety '
                                               'found'},
                                 {   'smiles': 'C1CCC(C1)CC#CC2=CC=C(C=C2)[C@@H]3[C@@H]4CN(CC(=O)N4[C@@H]3CO)C(=O)C5=CC=C(C=C5)F',
                                     'name': '(6R,7R,8S)-7-[4-(3-cyclopentylprop-1-ynyl)phenyl]-4-[(4-fluorophenyl)-oxomethyl]-8-(hydroxymethyl)-1,4-diazabicyclo[4.2.0]octan-2-one',
                                     'reason': 'No beta-D-glucose moiety '
                                               'found'},
                                 {   'smiles': 'S(=O)(C(SSCCC)CC)CCC',
                                     'name': 'Propyl 1-(propylsulfinyl)propyl '
                                             'disulfide',
                                     'reason': 'No beta-D-glucose moiety '
                                               'found'},
                                 {   'smiles': 'ClC=1C(=O)[C@@H]([C@@](O)(C/C=C\\CCCCC)C1)C[C@H](OC(=O)C)[C@@H](OC(=O)C)CCCC(OC)=O',
                                     'name': 'punaglandin 6',
                                     'reason': 'No beta-D-glucose moiety '
                                               'found'},
                                 {   'smiles': 'CCCCCCCCCCCCCCCCCCCCCCCC[C@H](O)C(=O)N[C@@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)[C@H](O)[C@@H](O)CCCCCCCCCCCCCC',
                                     'name': 'Ins-1-P-Cer(t18:0/2-OH-26:0)',
                                     'reason': 'No beta-D-glucose moiety '
                                               'found'},
                                 {   'smiles': 'CCC(C)(C)C(=O)OC1CC(C=C2C1[C@H]([C@H](C=C2)C)CC[C@@H]3CC(CC(=O)O3)O)C',
                                     'name': '2,2-dimethylbutanoic acid '
                                             '[(7S,8S)-8-[2-[(2R)-4-hydroxy-6-oxo-2-oxanyl]ethyl]-3,7-dimethyl-1,2,3,7,8,8a-hexahydronaphthalen-1-yl] '
                                             'ester',
                                     'reason': 'No beta-D-glucose moiety '
                                               'found'},
                                 {   'smiles': '[C@@H]1(C2=CC=CC=C2)[C@@H](N)C1',
                                     'name': '(1S,2R)-tranylcypromine',
                                     'reason': 'No beta-D-glucose moiety '
                                               'found'},
                                 {   'smiles': 'OC(=O)/C=C/C#CC#CC#N',
                                     'name': '(2E)-7-Cyanohept-2-en-4,6-diynoic '
                                             'acid',
                                     'reason': 'No beta-D-glucose moiety '
                                               'found'}],
    'sample_false_negatives': [   {   'smiles': 'CC[C@H]1O[C@@H](OCC\\C=C/C[C@@H]2[C@H](CCC2=O)CC(O)=O)[C@H](O)[C@@H](O)[C@@H]1O',
                                      'name': '12-hydroxyjasmonic acid '
                                              '12-O-beta-D-glucoside',
                                      'reason': 'No beta-D-glucose moiety '
                                                'found'},
                                  {   'smiles': 'O1[C@@H](OC[C@@H](O)[C@H](O)[C@H](O)[C@@H](NC(=O)C)CO)[C@H](NC(=O)C)[C@@H](O)[C@H](O)[C@H]1CO',
                                      'name': 'beta-D-GlcpNAc-(1->6)-GalNAc-ol',
                                      'reason': 'No beta-D-glucose moiety '
                                                'found'},
                                  {   'smiles': 'C[C@@H]1O[C@@H](Oc2cc(OC(C)=O)c3c(c2)oc(-c2ccc(OC(C)=O)c(OC(C)=O)c2)c(O[C@@H]2O[C@H](CC(C)=O)[C@@H](OC(C)=O)[C@H](OC(C)=O)[C@H]2OC(C)=O)c3=O)[C@H](OC(C)=O)[C@H](OC(C)=O)[C@H]1OC(C)=O',
                                      'name': 'petiolaroside decaacetate',
                                      'reason': 'No beta-D-glucose moiety '
                                                'found'},
                                  {   'smiles': '[H][C@@]1(O[C@@H](Oc2c(oc3cc(O)cc(O)c3c2=O)-c2ccc(O)c(O)c2)[C@H](O)[C@H]1O)[C@H](O)CO',
                                      'name': 'quercetin '
                                              '3-O-beta-D-glucofuranoside',
                                      'reason': 'No beta-D-glucose moiety '
                                                'found'}],
    'attempt': 0,
    'reasoning': 'To classify a molecule as a beta-D-glucoside, we need to '
                 'identify the presence of a beta-D-glucose moiety attached to '
                 'another molecule via a glycosidic bond. The key features to '
                 'look for are:\n'
                 '1. The presence of a glucose ring structure.\n'
                 '2. The glucose ring must be in the beta configuration (i.e., '
                 'the anomeric hydroxyl group is in the axial position).\n'
                 '3. The glucose must be attached to another molecule via a '
                 'glycosidic bond.\n'
                 '\n'
                 'Here is the Python program that implements this '
                 'classification:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 540,
    'num_false_positives': 100,
    'num_true_negatives': 810,
    'num_false_negatives': 4,
    'num_negatives': None,
    'precision': 0.84375,
    'recall': 0.9926470588235294,
    'f1': 0.9121621621621622,
    'accuracy': 0.9284731774415406,
    'negative_predictive_value': 0.995085995085995}