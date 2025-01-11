"""
Classifies: CHEBI:83970 cardiac glycoside
"""
"""
Classifies: cardiac glycoside
Steroid lactones containing sugar residues that act on the contractile force of the cardiac muscles.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cardiac_glycoside(smiles: str):
    """
    Determines if a molecule is a cardiac glycoside based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a cardiac glycoside, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check for steroid core (tetracyclic ring system)
    steroid_pattern = Chem.MolFromSmarts("[C]1[C][C]2[C]([C]1)[C]1[C][C][C]3[C]([C]1[C][C]2)[C][C][C]2[C][C][C][C][C]23")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid core found"
        
    # Check for butenolide/furanone ring
    lactone_pattern = Chem.MolFromSmarts("C1=CC(=O)OC1")
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone ring found"
        
    # Check for sugar moiety (pyranose ring with multiple OH groups)
    sugar_pattern = Chem.MolFromSmarts("[OX2H1][CX4H1][CX4H1][CX4H1][CX4H1][CX4H1]O")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar moiety found"
        
    # Check for glycosidic linkage (C-O-C between sugar and steroid)
    glycosidic_pattern = Chem.MolFromSmarts("[CX4][OX2][CX4]")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No glycosidic linkage found"
        
    # Count oxygen atoms (should have multiple due to OH groups and sugar)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 6:
        return False, "Too few oxygen atoms for cardiac glycoside"
        
    # Check molecular weight (typically >500 Da for cardiac glycosides)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for cardiac glycoside"
        
    # Count rings (should have multiple rings including steroid core + sugar)
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 5:  # steroid (4) + at least 1 sugar ring
        return False, "Too few rings for cardiac glycoside"
        
    # Count hydroxyl groups (should have multiple)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if hydroxyl_matches < 3:
        return False, "Too few hydroxyl groups"

    return True, "Contains steroid core with lactone ring and sugar residue(s)"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:83970',
                          'name': 'cardiac glycoside',
                          'definition': 'Steroid lactones containing sugar '
                                        'residues that act on the contractile '
                                        'force of the cardiac muscles.',
                          'parents': ['CHEBI:24400', 'CHEBI:26766'],
                          'xrefs': ['Wikipedia:Cardiac_glycoside'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'O1[C@H](CC2=C(C=3C=4C(C(O)=C5C3C[C@H](C)OC5)=C(OC)C=C(C4)OC)C6=CC(OC)=CC(=C6C(=C2C1)O)OC)C',
                                     'name': '(3S)-5-[(3S)-10-hydroxy-7,9-dimethoxy-3-methyl-3,4-dihydro-1H-benzo[g]isochromen-5-yl]-7,9-dimethoxy-3-methyl-3,4-dihydro-1H-benzo[g]isochromen-10-ol',
                                     'reason': 'No steroid core found'},
                                 {   'smiles': 'C[C@H]1CN([C@@H](COC2=C(C=CC(=C2)NC(=O)C)C(=O)N(C[C@@H]1OC)C)C)C(=O)C3CCC3',
                                     'name': 'N-[(5R,6S,9R)-8-[cyclobutyl(oxo)methyl]-5-methoxy-3,6,9-trimethyl-2-oxo-11-oxa-3,8-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]acetamide',
                                     'reason': 'No steroid core found'},
                                 {   'smiles': '[C@@]12(CCCC[C@@]1(C=C[C@@H]([C@@H]2CC[C@H](C[C@H](CC(=O)[O-])O)O)C)[H])[H]',
                                     'name': '4a,5-dihydro-ML-236C carboxylate',
                                     'reason': 'No steroid core found'},
                                 {   'smiles': 'S(=O)(CC1OC(N2C3=NC=NC(N)=C3N=C2)C(O)C1O)C',
                                     'name': "(S)-5'-Deoxy-5'-(methylsulfinyl)adenosine",
                                     'reason': 'No steroid core found'},
                                 {   'smiles': 'OC=1C(=C(C=2C=3C(NC2)=CC=CC3)C(=O)C(=O)C1C=4C=5C(NC4)=CC=CC5)CC=C(C)C',
                                     'name': 'Ochrindole D',
                                     'reason': 'No steroid core found'},
                                 {   'smiles': 'O(C1C(O)C(OC1OC2=C(OC=3C(C2=O)=C(O)C=C(OC4OC(C(O)C(O)C4O)C)C3)C5=CC=C(O)C=C5)CO)C6OCC(O)(C6O)CO',
                                     'name': 'Kaempferol '
                                             '3-apiosyl-(1->2)-alpha-L-arabinofuranoside-7-rhamnoside',
                                     'reason': 'No steroid core found'},
                                 {   'smiles': '[O-]C(=O)OON=O',
                                     'name': 'nitrosoperoxycarbonate(1-)',
                                     'reason': 'No steroid core found'},
                                 {   'smiles': 'P(OCC(O)COC(=O)CCCCCCCCCCCCCCCCC)(O)(O)=O',
                                     'name': 'LysoPA(18:0/0:0)',
                                     'reason': 'No steroid core found'},
                                 {   'smiles': 'O=C1C=2C(=O)OC34C2OC(C(O)C3C=CC(C(CC1C)C)=O)C(C)C4',
                                     'name': 'Atrop-Abybetaomicin C',
                                     'reason': 'No steroid core found'},
                                 {   'smiles': 'P(OC[C@H](OC(=O)CCCCCCCCCCCCCCCCCC)CO/C=C\\CCCCCCCCCCCCCCCC)(OC[C@@H](O)CO)(O)=O',
                                     'name': 'PG(P-18:0/19:0)',
                                     'reason': 'No steroid core found'}],
    'sample_false_negatives': [   {   'smiles': 'O[C@@]12[C@]3(C([C@@]4(C(CC3)C[C@]5(O[C@@]6(O[C@@H](C[C@@H](O)[C@@]6(O[C@@]5(C4)[H])O)C)[H])[H])C=O)CC[C@@]1([C@H](CC2)C=7COC(=O)C7)C)[H]',
                                      'name': '(3R,5S,7R,9R,10S,12R,14R,18R,19R,22S,23R)-9,10,22-trihydroxy-7,18-dimethyl-19-(5-oxo-2H-furan-3-yl)-4,6,11-trioxahexacyclo[12.11.0.03,12.05,10.015,23.018,22]pentacosane-14-carbaldehyde',
                                      'reason': 'No steroid core found'},
                                  {   'smiles': 'C[C@@H]1O[C@@H](O[C@H]2CC[C@@]3(CO)[C@@H](CC[C@@H]4[C@@H]3[C@@H](O)C[C@]3(C)[C@H](CC[C@]43O)C3=CC(=O)OC3)C2)[C@H](O)[C@H](O)[C@H]1O',
                                      'name': 'Panogenin-3-o-alpha-L-rhamnopyranoside',
                                      'reason': 'No steroid core found'},
                                  {   'smiles': 'CO[C@@H]1[C@@H](O)[C@H](C)O[C@@H](O[C@H]2CC[C@@]3(C=O)[C@H](CC[C@@H]4[C@@H]3[C@@H]3O[C@@H]3[C@]3(C)[C@H](CC[C@]43O)C3=CC(=O)OC3)C2)[C@H]1OC(C)=O',
                                      'name': 'Cerbertin',
                                      'reason': 'No steroid core found'},
                                  {   'smiles': 'O[C@@]12[C@]3([C@@]([C@@]4([C@](O)(CC3)C[C@@H](O[C@@H]5O[C@@H]([C@@H](O[C@@H]6O[C@@H]([C@@H](O)[C@H](O)[C@H]6O)CO)[C@@H](O)C5)C)CC4)C)([C@H](O)C[C@@]1([C@H](CC2)C=7COC(=O)C7)C)[H])[H]',
                                      'name': 'Corchorusoside A',
                                      'reason': 'No steroid core found'},
                                  {   'smiles': '[H][C@]12CC[C@]3([H])[C@]([H])(C[C@@H](OC(C)=O)[C@]4(C)[C@H](CC[C@]34O)C3=CC(=O)OC3)[C@@]1(C)CC[C@@H](C2)O[C@H]1C[C@H](O)[C@H](O[C@H]2C[C@H](O)[C@H](O[C@H]3C[C@H](O)[C@H](O)[C@@H](C)O3)[C@@H](C)O2)[C@@H](C)O1',
                                      'name': '12-O-acetyldigoxin',
                                      'reason': 'No steroid core found'},
                                  {   'smiles': 'OC12C3C(C4(C(O)(CC3)CC(OC5OC(C(O)C(O)C5O)C)CC4)C)C(O)CC1(C(CC2)C=6COC(=O)C6)C',
                                      'name': 'Bipindogulomethyloside',
                                      'reason': 'No steroid core found'},
                                  {   'smiles': 'CO[C@H]1[C@@H](O)[C@@H](C)O[C@@H](O[C@H]2CC[C@@]3(C)[C@H](CC[C@@H]4[C@@H]3C(=O)[C@@H](O)[C@]3(C)[C@H](CC[C@]43O)C3=CC(=O)OC3)C2)[C@@H]1O',
                                      'name': 'Musaroside',
                                      'reason': 'No steroid core found'},
                                  {   'smiles': 'C[C@@H]1C[C@@H](O)[C@]2(O)O[C@@H]3C[C@@]4(C=O)[C@@H](CC[C@@H]5[C@@H]4CC[C@]4(C)[C@H](CC[C@]54O)C4=CC(=O)OC4)C[C@H]3O[C@@H]2O1',
                                      'name': 'Calactin',
                                      'reason': 'No steroid core found'},
                                  {   'smiles': 'O[C@@]12[C@]3([C@@]([C@@]4([C@](O)(CC3)C[C@@H](O[C@@H]5O[C@@H]([C@@H](O[C@@H]6O[C@@H]([C@@H](O)[C@H](O)[C@H]6O)CO)[C@@H](O)C5)C)CC4)C)(CC[C@@]1([C@H](CC2)C=7COC(=O)C7)C)[H])[H]',
                                      'name': 'Corchorusoside C',
                                      'reason': 'No steroid core found'},
                                  {   'smiles': 'CO[C@H]1C[C@H](O[C@H]2CC[C@@]3(C)[C@H](CC[C@@H]4[C@@H]3[C@H](O)C[C@]3(C)[C@H](CC[C@]43O)C3=CC(=O)OC3)C2)O[C@@H](C)[C@H]1O',
                                      'name': 'Divostroside',
                                      'reason': 'No steroid core found'}],
    'attempt': 0,
    'reasoning': 'Let me help create a program to classify cardiac glycosides. '
                 'From the definition and examples, I can identify several key '
                 'structural features:\n'
                 '\n'
                 '1. A steroid core (tetracyclic ring system)\n'
                 '2. A lactone ring (specifically a butenolide/furanone)\n'
                 '3. At least one sugar residue attached via glycosidic bond\n'
                 '4. Hydroxyl groups are common\n'
                 '5. Often has an unsaturated lactone ring at C-17\n'
                 '\n'
                 "Here's the code to identify cardiac glycosides:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 142214,
    'num_false_negatives': 86,
    'num_negatives': None,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': 0.9993956430077301,
    'negative_predictive_value': 0.9993956430077301}