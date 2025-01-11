"""
Classifies: CHEBI:17002 cholesteryl ester
"""
"""
Classifies: cholesteryl ester
A sterol ester obtained by formal condensation of the carboxy group of any 
carboxylic acid with the 3-hydroxy group of cholesterol.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cholesteryl_ester(smiles: str):
    """
    Determines if a molecule is a cholesteryl ester based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a cholesteryl ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Look for steroid core (four fused rings)
    steroid_core = Chem.MolFromSmarts("[#6]12[#6][#6][#6]3[#6][#6][#6]4[#6][#6][#6][#6]([#6][#6][#6]2[#6][#6][#6]1)[#6]34")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Look for exactly one ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Check for cholesterol-specific features
    # Including the double bond at C5-C6 and methyl groups
    cholesterol_pattern = Chem.MolFromSmarts("[#6]-[#6]1-[#6]-[#6]-[#6]2=[#6]-[#6]-[#6]-[#6]3-[#6]-[#6]-[#6]-[#6]-[#6]-3-[#6]-[#6]-2-[#6]-1")
    if not mol.HasSubstructMatch(cholesterol_pattern):
        return False, "Missing characteristic cholesterol features"

    # Verify ester is attached at position 3 (the only oxygen should be part of the ester)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count != 2:  # One for ester C=O and one for ester C-O-
        return False, f"Found {oxygen_count} oxygens, should be exactly 2 for cholesteryl ester"

    # Check molecular properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:  # Cholesteryl esters are typically >500 Da
        return False, "Molecular weight too low for cholesteryl ester"

    # Count carbons - cholesteryl esters typically have >30 carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 29:
        return False, "Too few carbons for cholesteryl ester"

    # Check for reasonable chain length of fatty acid part
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Acyl chain too short for typical cholesteryl ester"

    return True, "Contains cholesterol core with ester-linked fatty acid at position 3"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:17002',
                          'name': 'cholesteryl ester',
                          'definition': 'A sterol ester obtained by formal '
                                        'condensation of the carboxy group of '
                                        'any carboxylic acid with the '
                                        '3-hydroxy group of cholesterol.',
                          'parents': ['CHEBI:35915'],
                          'xrefs': ['KEGG:C02530', 'PMID:2474544'],
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
    'sample_false_negatives': [   {   'smiles': 'O([C@@H]1CC=2[C@@](C3C(C4[C@@](C(CC4)[C@@H](CCCC(C)C)C)(CC3)C)CC2)(CC1)C)C(CCCCCCCCCCC/C=C\\CCCCCCCC)=O',
                                      'name': 'CE(22:1(13Z))',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': 'O([C@@H]1CC=2[C@@]([C@@]3([C@]([C@]4([C@@]([C@](CC4)([C@@H](CCCC(C)C)C)[H])(CC3)C)[H])(CC2)[H])[H])(CC1)C)C(=O)CCCCCCCCCCCCCCCCCCC',
                                      'name': 'cholest-5-en-3beta-yl '
                                              'eicosanoate',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': 'CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCC(=O)O[C@H]1CC[C@]2(C)[C@H]3CC[C@]4(C)[C@H](CC[C@H]4[C@@H]3CC=C2C1)[C@H](C)CCCC(C)C',
                                      'name': 'cholesteryl '
                                              '(7Z,10Z,13Z,16Z,19Z)-docosapentaenoate',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': 'O([C@@H]1CC=2[C@@]([C@@]3([C@]([C@]4([C@@]([C@](CC4)([C@@H](CCCC(C)C)C)[H])(CC3)C)[H])(CC2)[H])[H])(CC1)C)C(=O)CCCCCCCCCCC',
                                      'name': 'Cholesteryl laurate',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': 'CCCCCCCCCCCCCCC(=O)O[C@H]1CC[C@]2(C)[C@H]3CC[C@]4(C)[C@H](CC[C@H]4[C@@H]3CC=C2C1)[C@H](C)CCCC(C)C',
                                      'name': 'cholesteryl pentadecanoate',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': 'O([C@@H]1CC=2[C@@](C3C(C4[C@@](C(CC4)[C@@H](CCCC(C)C)C)(CC3)C)CC2)(CC1)C)C(CCCCCCCCC/C=C\\C/C=C\\CCCCC)=O',
                                      'name': 'CE(20:2(6Z,9Z))',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': 'O([C@@H]1CC=2[C@@]([C@@]3([C@]([C@]4([C@@]([C@](CC4)([C@@H](CCCC(C)C)C)[H])(CC3)C)[H])(CC2)[H])[H])(CC1)C)C(=O)CCCCCCCCCCC/C=C\\CCCCCCCC',
                                      'name': 'cholest-5-en-3b-yl '
                                              '(13Z-docosenoate)',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': '[H][C@@]1(CC[C@@]2([H])[C@]3([H])CC=C4C[C@H](CC[C@]4(C)[C@@]3([H])CC[C@]12C)OC(=O)CCCCCCCCC\\C=C/CCCCCC)[C@H](C)CCCC(C)C',
                                      'name': '(11Z-octadecenoyl)-cholesterol',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': 'C1[C@@H](CC=2[C@](C1)([C@@]3([C@@](CC2)([C@]4([C@@](CC3)(C)C(CC4)[C@@H](CCCC(C)C)C)[H])[H])[H])C)OC(=O)CCCCCCC/C=C/C/C=C/CCC',
                                      'name': 'CE(16:2)',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': 'CC\\C=C/C\\C=C/C\\C=C/CCCCCCCC(=O)O[C@H]1CC[C@]2(C)[C@H]3CC[C@]4(C)[C@H](CC[C@H]4[C@@H]3CC=C2C1)[C@H](C)CCCC(C)C',
                                      'name': 'cholesteryl linolenate',
                                      'reason': 'No steroid core structure '
                                                'found'}],
    'attempt': 0,
    'reasoning': 'Let me help create a program to classify cholesteryl esters. '
                 'These compounds consist of a cholesterol backbone with a '
                 'fatty acid attached via an ester bond at the 3-hydroxy '
                 'position.\n'
                 '\n'
                 "Here's my step-by-step approach:\n"
                 '1. Check for the cholesterol core structure\n'
                 '2. Verify the presence of exactly one ester group\n'
                 '3. Verify the ester is attached at the correct position\n'
                 '4. Check for reasonable molecular properties',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 142259,
    'num_false_negatives': 41,
    'num_negatives': None,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': 0.9997118763176388,
    'negative_predictive_value': 0.9997118763176388}