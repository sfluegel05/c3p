"""
Classifies: CHEBI:35436 D-glucoside
"""
"""
Classifies: CHEBI:16238 D-glucoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_D_glucoside(smiles: str):
    """
    Determines if a molecule is a D-glucoside based on its SMILES string.
    A D-glucoside is a molecule where a D-glucose moiety is attached to another group via a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a D-glucoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the D-glucose pattern (beta-D-glucopyranose)
    glucose_pattern = Chem.MolFromSmarts("[C@H]1([C@H]([C@@H]([C@H]([C@@H](O1)CO)O)O)O)")
    if not mol.HasSubstructMatch(glucose_pattern):
        return False, "No D-glucose moiety found"

    # Check for a glycosidic bond (oxygen connected to the anomeric carbon of glucose)
    glycosidic_bond_pattern = Chem.MolFromSmarts("[C@H]1([C@H]([C@@H]([C@H]([C@@H](O1)CO)O)O)O)-O")
    if not mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "No glycosidic bond found"

    # Check that the glucose moiety is attached to another group
    # We do this by ensuring that the glucose moiety is not free (i.e., it has at least one bond to a non-hydrogen atom)
    glucose_atoms = mol.GetSubstructMatch(glucose_pattern)
    glucose_mol = Chem.PathToSubmol(mol, glucose_atoms)
    for atom in glucose_mol.GetAtoms():
        if atom.GetAtomicNum() == 8:  # Oxygen atom
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() != 1:  # Not hydrogen
                    return True, "Contains a D-glucose moiety attached via a glycosidic bond"

    return False, "Glucose moiety is not attached to another group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35436',
                          'name': 'D-glucoside',
                          'definition': 'Any glucoside in which the glycoside '
                                        'group is derived from D-glucose.',
                          'parents': ['CHEBI:24278'],
                          'xrefs': ['KEGG:C01798'],
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
               'CCCCCCCCCCCCCCCCC(O)C(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCC(C)C '
               'NAME: '
               'N-(2-hydroxyoctadecanoyl)-1-O-beta-D-glucosyl-15-methylhexadecasphing-4-enine '
               'REASON: MISSED Glucose moiety is not attached to another '
               'group\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCCCCCCC(O)C(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)[C@H](O)CCCCCCCCCCC(C)C '
               'NAME: '
               'N-(2-hydroxytetracosanoyl)-1-O-beta-D-glucosyl-4-hydroxy-15-methylhexadecasphinganine '
               'REASON: MISSED Glucose moiety is not attached to another '
               'group\n'
               ' * SMILES: '
               'O([C@@H]1O[C@@H]([C@H]([C@@H]([C@H]1O)O)O)COC(CC(=O)O)=O)C=2C=CC3=C(C2)C(=CN3)C(=O)O '
               'NAME: 5-(6-O-malonyl-beta-D-glucosyloxy)-indole-3-carboxylic '
               'acid REASON: MISSED Glucose moiety is not attached to another '
               'group\n'
               ' * SMILES: '
               'O(C=1C=CC(=CC1OC)[C@H]2OC3=CC(=CC(=C3C(C2)=O)O)O[C@H]4[C@@H]([C@H]([C@@H]([C@H](O4)CO)O)O)O)[C@@H]5OC[C@@]([C@@]5([H])O)(CO)O '
               'NAME: viscumneoside I REASON: MISSED Glucose moiety is not '
               'attached to another group\n'
               ' * SMILES: '
               'OC[C@H]1O[C@@H](Oc2ccc(COC(=O)c3ccc(O)c(O)c3)cc2)[C@H](O)[C@@H](O)[C@@H]1O '
               'NAME: amburoside A REASON: MISSED Glucose moiety is not '
               'attached to another group\n'
               ' * SMILES: '
               'OC[C@H]1O[C@@H](OC(=O)c2cc(O)c(O)c(O)c2)[C@H](O)[C@@H](O)[C@@H]1O '
               'NAME: 1-O-galloyl-beta-D-glucose REASON: MISSED Glucose moiety '
               'is not attached to another group\n'
               ' * SMILES: '
               'COc1cc(ccc1O)-c1[o+]c2cc(O)cc(O)c2cc1O[C@@H]1O[C@H](COC(C)=O)[C@@H](O)[C@H](O)[C@H]1O '
               'NAME: peonidin 3-O-(6-O-acetyl-beta-D-glucoside) REASON: '
               'MISSED Glucose moiety is not attached to another group\n'
               ' * SMILES: '
               '[H][C@]1(OC[C@](O)(COC(=O)c2cc(OC)c(O)c(OC)c2)[C@H]1O)O[C@@H]1[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]1Oc1cc(OC)c(OC)c(OC)c1 '
               'NAME: albibrissinoside A REASON: MISSED Glucose moiety is not '
               'attached to another group\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)[C@H](O)CCCCCCCCCCC(C)C '
               'NAME: '
               'N-henicosanoyl-1-O-beta-D-glucosyl-4-hydroxy-15-methylhexadecasphinganine '
               'REASON: MISSED Glucose moiety is not attached to another '
               'group\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)[C@H](O)CCCCCCCCCCC(C)C '
               'NAME: '
               'N-pentacosanoyl-1-O-beta-D-glucosyl-4-hydroxy-15-methylhexadecasphinganine '
               'REASON: MISSED Glucose moiety is not attached to another '
               'group\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': '[Si](OC(=O)CCCCCCCCCCC)(C)(C)C',
                                     'name': 'Dodecanoic acid, trimethylsilyl '
                                             'ester',
                                     'reason': 'No D-glucose moiety found'},
                                 {   'smiles': 'C[C@H]1CN(C(=O)C2=CC=CC=C2C3=CC=CC=C3CO[C@H]1CN(C)C(=O)NC4=C(C=C(C=C4)OC)OC)[C@H](C)CO',
                                     'name': 'LSM-30778',
                                     'reason': 'No D-glucose moiety found'},
                                 {   'smiles': 'C(C([C@H]([C@@H]([C@H]([C@@H](CO)O)O)O)O)=O)O',
                                     'name': 'D-ido-heptulose',
                                     'reason': 'No D-glucose moiety found'},
                                 {   'smiles': 'O=C1C=C(CC(O)C)O[C@]1(/C=C/C=C/CC)C',
                                     'name': 'Terrefuranone',
                                     'reason': 'No D-glucose moiety found'},
                                 {   'smiles': 'C1CC1CNC(=O)C[C@@H]2C[C@@H]3[C@H]([C@H](O2)CO)OC4=C3C=C(C=C4)NC(=O)CC5=CC=NC=C5',
                                     'name': '2-[(1R,3S,4aS,9aR)-1-(hydroxymethyl)-6-[(1-oxo-2-pyridin-4-ylethyl)amino]-3,4,4a,9a-tetrahydro-1H-pyrano[3,4-b]benzofuran-3-yl]-N-(cyclopropylmethyl)acetamide',
                                     'reason': 'No D-glucose moiety found'},
                                 {   'smiles': 'COC1=CC=C(C=C1)C2=NN3C(=C2)N=C4CCCC4=C3N5CCC6(CC5)OCCO6',
                                     'name': 'LSM-28273',
                                     'reason': 'No D-glucose moiety found'},
                                 {   'smiles': '[O-]S(C[C@H]([C@@H](O)C=O)O)(=O)=O',
                                     'name': '4-deoxy-4-sulfo-D-erythrose(1-)',
                                     'reason': 'No D-glucose moiety found'},
                                 {   'smiles': 'O[C@@H](CCCCC)\\C=C\\CCCCCCCCCC(O)=O',
                                     'name': '13S-HOME(11E)',
                                     'reason': 'No D-glucose moiety found'},
                                 {   'smiles': 'CC(CCCCCCC/C=C\\C(=O)NO)C',
                                     'name': '(Z)-N-hydroxy-11-methyldodec-2-enamide',
                                     'reason': 'No D-glucose moiety found'},
                                 {   'smiles': 'CCCC1=CC(=CS1)C2=NNC(=S)N2C(C)C3=CC=CC=C3',
                                     'name': '4-(1-phenylethyl)-3-(5-propyl-3-thiophenyl)-1H-1,2,4-triazole-5-thione',
                                     'reason': 'No D-glucose moiety found'}],
    'sample_false_negatives': [   {   'smiles': 'CCN1C(=O)CC(SC[C@H](NC(C)=O)C(=O)N[C@@H]2[C@@H](O)[C@H](O)[C@@H](CO)O[C@@H]2O[C@H]2[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]2O)C1=O',
                                      'name': 'N-ethylsuccinimide-mycothiol '
                                              'conjugate',
                                      'reason': 'No D-glucose moiety found'},
                                  {   'smiles': 'CC[C@H]1O[C@@H](OCC\\C=C/C[C@@H]2[C@H](CCC2=O)CC(O)=O)[C@H](O)[C@@H](O)[C@@H]1O',
                                      'name': '12-hydroxyjasmonic acid '
                                              '12-O-beta-D-glucoside',
                                      'reason': 'No D-glucose moiety found'},
                                  {   'smiles': 'NC[C@H]1O[C@H](O[C@@H]2[C@@H](N)C[C@@H](N)[C@H](O)[C@H]2O)[C@H](N)[C@@H](O)[C@@H]1O',
                                      'name': 'neamine',
                                      'reason': 'No D-glucose moiety found'},
                                  {   'smiles': '[H]C(=O)SC[C@H](NC(C)=O)C(=O)N[C@@H]1[C@@H](O)[C@H](O)[C@@H](CO)O[C@@H]1O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O',
                                      'name': 'S-formylmycothiol',
                                      'reason': 'No D-glucose moiety found'},
                                  {   'smiles': 'N[C@@H]1[C@@H](O)[C@H](O)[C@@H](CO)O[C@@H]1O[C@@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@@H](O)[C@H]1OP(O)(=O)OCCCCCCS',
                                      'name': '6-O-(2-amino-2-deoxy-alpha-D-glucosyl)-1D-myo-inositol '
                                              '1-(6-mercaptohexyl)phosphate',
                                      'reason': 'No D-glucose moiety found'},
                                  {   'smiles': 'O1[C@@H](OC[C@@H](O)[C@H](O)[C@H](O)[C@@H](NC(=O)C)CO)[C@H](NC(=O)C)[C@@H](O)[C@H](O)[C@H]1CO',
                                      'name': 'beta-D-GlcpNAc-(1->6)-GalNAc-ol',
                                      'reason': 'No D-glucose moiety found'},
                                  {   'smiles': 'C[C@@H]1O[C@@H](Oc2cc(OC(C)=O)c3c(c2)oc(-c2ccc(OC(C)=O)c(OC(C)=O)c2)c(O[C@@H]2O[C@H](CC(C)=O)[C@@H](OC(C)=O)[C@H](OC(C)=O)[C@H]2OC(C)=O)c3=O)[C@H](OC(C)=O)[C@H](OC(C)=O)[C@H]1OC(C)=O',
                                      'name': 'petiolaroside decaacetate',
                                      'reason': 'No D-glucose moiety found'},
                                  {   'smiles': 'CC(=O)N[C@@H]1[C@@H](O)[C@H](O)[C@@H](CO)O[C@@H]1O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](OP(O)(O)=O)[C@H]1O',
                                      'name': '1D-myo-inositol '
                                              '2-acetamido-2-deoxy-alpha-D-glucopyranoside '
                                              '3-phosphate',
                                      'reason': 'No D-glucose moiety found'},
                                  {   'smiles': '[H][C@@]1(O[C@@H](Oc2c(oc3cc(O)cc(O)c3c2=O)-c2ccc(O)c(O)c2)[C@H](O)[C@H]1O)[C@H](O)CO',
                                      'name': 'quercetin '
                                              '3-O-beta-D-glucofuranoside',
                                      'reason': 'No D-glucose moiety found'},
                                  {   'smiles': 'C([C@H](CS)N)(N[C@@H]1[C@H]([C@@H]([C@H](O[C@@H]1O[C@@H]2[C@@H]([C@@H]([C@H]([C@@H]([C@H]2O)O)O)O)O)CO)O)O)=O',
                                      'name': '1D-myo-inositol '
                                              '2-(L-cysteinylamino)-2-deoxy-alpha-D-glucopyranoside',
                                      'reason': 'No D-glucose moiety found'}],
    'attempt': 1,
    'reasoning': 'The previous code failed to correctly classify D-glucosides '
                 'because it incorrectly assumed that the glucose moiety must '
                 'be attached to another group via a glycosidic bond. However, '
                 'the error messages indicate that the glucose moiety was not '
                 'recognized as being attached to another group, even when it '
                 'was. This suggests that the pattern matching for the '
                 'glycosidic bond and the free glucose pattern was too '
                 'restrictive.\n'
                 '\n'
                 'To improve the classification, we need to:\n'
                 '1. Ensure that the glucose moiety is present.\n'
                 '2. Check that the glucose moiety is attached to another '
                 'group via a glycosidic bond.\n'
                 '3. Avoid false negatives by not requiring the glucose moiety '
                 'to be free.\n'
                 '\n'
                 'Here’s the revised code:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 583,
    'num_false_positives': 100,
    'num_true_negatives': 787,
    'num_false_negatives': 16,
    'num_negatives': None,
    'precision': 0.8535871156661786,
    'recall': 0.9732888146911519,
    'f1': 0.9095163806552262,
    'accuracy': 0.9219380888290714,
    'negative_predictive_value': 0.9800747198007472}