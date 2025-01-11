"""
Classifies: CHEBI:1624 3-oxo-5beta-steroid
"""
"""
Classifies: CHEBI:35341 3-oxo-5beta-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_5beta_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5beta-steroid based on its SMILES string.
    A 3-oxo-5beta-steroid has a ketone at position 3 and beta configuration at position 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-5beta-steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for basic steroid core (4 fused rings)
    steroid_core = Chem.MolFromSmarts("[C]1[C][C]2[C][C][C]3[C][C][C]4[C][C][C][C]4[C]3[C][C]2[C]1")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Check for ketone at position 3
    # The pattern looks for a carbonyl (C=O) connected to two carbons in the first ring
    ketone_pattern = Chem.MolFromSmarts("[CH2][C](=O)[CH2]")
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "No ketone group at position 3"

    # Check for 5-beta configuration
    # In 5-beta steroids, the A/B ring junction is trans
    # We can look for the specific stereochemistry pattern
    # The pattern checks for the characteristic trans fusion at the A/B ring junction
    # with the hydrogen at position 5 being up (beta)
    ab_junction_pattern = Chem.MolFromSmarts("[C]1[CH2][C](=O)[CH2][C@@H]2[CH2]")
    if not mol.HasSubstructMatch(ab_junction_pattern):
        return False, "Not a 5-beta configuration"

    # Additional validation: Check carbon count (steroids typically have 19-27 carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 19 or carbon_count > 30:
        return False, f"Carbon count ({carbon_count}) outside typical steroid range (19-30)"

    # If we get here, it's a 3-oxo-5beta-steroid
    return True, "Molecule contains steroid core with 3-oxo group and 5-beta configuration"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:1624',
                          'name': '3-oxo-5beta-steroid',
                          'definition': 'Any 3-oxo steroid that has beta- '
                                        'configuration at position 5.',
                          'parents': ['CHEBI:136889', 'CHEBI:47788'],
                          'xrefs': [   'KEGG:C02797',
                                       'MetaCyc:3-Oxo-5-Beta-Steroids'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'S(=O)(=O)(C1=C2C(=CC3=C1NC=4C=CC=CC34)[C@@]5([C@H]([C@](C(=O)O)([C@@H](O)CC5)C)CC2)C)C6=CC7=C(NC8=C7C=C9[C@@]%10([C@H]([C@](C(=O)O)([C@@H](O)CC%10)C)CCC9=C8)C)C=C6',
                                     'name': 'Sulfadixiamycin C',
                                     'reason': 'No steroid core structure '
                                               'found'},
                                 {   'smiles': 'CNC(O)=O',
                                     'name': 'methylcarbamic acid',
                                     'reason': 'No steroid core structure '
                                               'found'},
                                 {   'smiles': 'CCNC(=O)NC1=CC2=C(C=C1)OC[C@H]3[C@@H](CC[C@H](O3)CC(=O)N[C@@H](C)C4=CC=CC=C4)N(C2=O)C',
                                     'name': '2-[(2S,4aR,12aR)-8-(ethylcarbamoylamino)-5-methyl-6-oxo-2,3,4,4a,12,12a-hexahydropyrano[2,3-c][1,5]benzoxazocin-2-yl]-N-[(1S)-1-phenylethyl]acetamide',
                                     'reason': 'No steroid core structure '
                                               'found'},
                                 {   'smiles': 'O([C@H]1O[C@@H]([C@@H](O[C@@H]2O[C@@H]([C@@H](O[C@@H]3O[C@@H]([C@H](O)[C@H](O)[C@H]3O)CO[C@]4(O[C@H]([C@H](NC(=O)C)[C@@H](O)C4)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]2NC(=O)C)CO)[C@H](O)[C@@H]1O[C@@H]5O[C@@H]([C@@H](O[C@@H]6O[C@@H]([C@H](O)[C@H](O)[C@H]6O)CO[C@]7(O[C@H]([C@H](NC(=O)C)[C@@H](O)C7)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]5NC(=O)C)CO)CO)[C@H]8[C@H](O)[C@H](O[C@@H](O[C@H]9[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]9CO)O[C@H]%10[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]%10CO[C@@H]%11O[C@H]([C@@H](O)[C@@H](O)[C@@H]%11O)C)O)[C@H]8O)CO[C@H]%12O[C@@H]([C@@H](O)[C@H](O)[C@@H]%12O[C@@H]%13O[C@@H]([C@@H](O[C@@H]%14O[C@@H]([C@H](O)[C@H](O[C@]%15(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%15)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]%14O)CO)[C@H](O)[C@H]%13NC(=O)C)CO)CO[C@@H]%16O[C@@H]([C@@H](O[C@@H]%17O[C@@H]([C@H](O)[C@H](O)[C@H]%17O)CO[C@]%18(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%18)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]%16NC(=O)C)CO',
                                     'name': 'CID 91851985',
                                     'reason': 'No steroid core structure '
                                               'found'},
                                 {   'smiles': 'O(C(=O)C(C1C(CN2C(C1)C=3NC=4C(C3CC2)=CC=CC4)CC)=COC)C',
                                     'name': 'Methyl '
                                             '2-(3-ethyl-1,2,3,4,6,7,12,12b-octahydroindolo[2,3-a]quinolizin-2-yl)-3-methoxyprop-2-enoate',
                                     'reason': 'No steroid core structure '
                                               'found'},
                                 {   'smiles': 'O[C@H](/C=C/C=C/C=C/[C@H](O)[C@H](O)C=C)[C@H](O)/C=C/C',
                                     'name': 'Separacene C',
                                     'reason': 'No steroid core structure '
                                               'found'},
                                 {   'smiles': 'C[C@@H]1O[C@@H](O[C@@H]2[C@@H](CO)O[C@@H](O[C@@H]3[C@@H](O)C(O)O[C@H](CO)[C@@H]3O)[C@H](NC(C)=O)[C@H]2O[C@@H]2O[C@H](CO)[C@H](O)[C@H](OS(O)(=O)=O)[C@H]2O)[C@@H](O)[C@H](O)[C@@H]1O',
                                     'name': 'alpha-L-Fucp-(1->4)-[beta-D-Galp3S-(1->3)]-beta-D-GlcpNAc-(1->3)-D-Galp',
                                     'reason': 'No steroid core structure '
                                               'found'},
                                 {   'smiles': 'C1=CC=CC2=C1C(N([C@H](C(N2)=O)CC=3C=CC(=CC3)OC)C)=O',
                                     'name': "(S)-4'-methoxycyclopeptine",
                                     'reason': 'No steroid core structure '
                                               'found'},
                                 {   'smiles': 'O=C(N[C@@H](CC=1C=2C(NC1)=CC=CC2)C(O)=O)[C@@H](NC(=O)[C@@H](N)C(C)C)C(C)C',
                                     'name': 'Val-Val-Trp',
                                     'reason': 'No steroid core structure '
                                               'found'},
                                 {   'smiles': 'C=1C(=C(C=CC1/C=C/CO)OC(CO)C(O)C=2C=C(C(=CC2)O)OC)OC',
                                     'name': 'guaiacylglycerol beta-coniferyl '
                                             'ether',
                                     'reason': 'No steroid core structure '
                                               'found'}],
    'sample_false_negatives': [   {   'smiles': 'C1[C@@]2([C@]3(C[C@@H]([C@]4([C@]([C@@]3(CC[C@@]2(CC(C1)=O)[H])[H])(CC[C@@]4([C@@H](CCC(O)=O)C)[H])[H])C)O)[H])C',
                                      'name': '12alpha-hydroxy-3-oxo-5beta-cholan-24-oic '
                                              'acid',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': '[H][C@]12CC[C@@]3([H])[C@]4([H])CC[C@](O)(C(=O)CO)[C@@]4(C)CC(=O)[C@]3([H])[C@@]1(C)CCC(=O)C2',
                                      'name': '17,21-dihydroxy-5beta-pregnane-3,11,20-trione',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': 'C[C@]12CC[C@H]3[C@@H](CC[C@@H]4CC(=O)CC[C@]34C)[C@@H]1CC[C@H]2O',
                                      'name': '5beta-dihydroepitestosterone',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': '[H][C@]12C[C@@H](O)[C@]3([H])[C@]([H])(CC[C@]4(C)[C@]([H])(CC[C@@]34[H])[C@H](C)CCCC(C)C)[C@@]1(C)CCC(=O)C2',
                                      'name': '7alpha-hydroxy-5beta-cholestan-3-one',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': 'C1C(C[C@@]2([C@](C1)([C@@]3([C@@](CC2)([C@@]4([H])[C@@](C[C@@H]3O)(C)[C@](CC4)(C(CO)=O)O)[H])[H])C)[H])=O',
                                      'name': '5beta-dihydrocortisol',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': 'C[C@]12CC[C@H]3[C@@H](CC[C@@H]4CC(=O)CC[C@]34C)[C@@H]1CC[C@@H]2C(=O)CO',
                                      'name': '5beta-dihydrodeoxycorticosterone',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': '[H][C@]12C[C@@H](O)[C@@]3([H])[C@]4([H])CC[C@]([H])([C@H](C)CCC(O)C(C)C)[C@@]4(C)CC[C@]3([H])[C@@]1(C)CCC(=O)C2',
                                      'name': '7alpha,24-dihydroxy-5beta-cholestan-3-one',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': 'C[C@]12C[C@H](O)[C@H]3[C@@H](CC[C@@H]4CC(=O)CC[C@]34C)[C@@H]1CC[C@@H]2C(=O)CO',
                                      'name': '5beta-dihydrocorticosterone',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': '[H][C@]12C[C@@H](O)[C@@]3([H])[C@]4([H])CC[C@]([H])([C@H](C)CC[C@H](O)C(C)C)[C@@]4(C)[C@@H](O)C[C@]3([H])[C@@]1(C)CCC(=O)C2',
                                      'name': '(24S)-7alpha,12alpha,24-trihydroxy-5beta-cholestan-3-one',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': 'CC(=O)[C@H]1CC[C@H]2[C@@H]3CC[C@@H]4CC(=O)CC[C@]4(C)[C@H]3[C@@H]3O[C@@H]3[C@]12C',
                                      'name': '11alpha,12alpha-Epoxy-5beta-pregnane-3,20-dione',
                                      'reason': 'No steroid core structure '
                                                'found'}],
    'attempt': 0,
    'reasoning': 'Let me help create a program to classify '
                 '3-oxo-5beta-steroids. From analyzing the examples and '
                 'definition, we need to check for:\n'
                 '\n'
                 '1. A steroid core structure (4 fused rings)\n'
                 '2. A ketone group (C=O) at position 3\n'
                 '3. Beta configuration at position 5 (H pointing up in the '
                 'traditional steroid representation)\n'
                 '\n'
                 "Here's the program:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 142273,
    'num_false_negatives': 27,
    'num_negatives': None,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': 0.9998102600140548,
    'negative_predictive_value': 0.9998102600140548}