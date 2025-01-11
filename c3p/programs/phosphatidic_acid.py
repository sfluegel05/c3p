"""
Classifies: CHEBI:16337 phosphatidic acid
"""
"""
Classifies: CHEBI:16337 phosphatidic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidic_acid(smiles: str):
    """
    Determines if a molecule is a phosphatidic acid based on its SMILES string.
    A phosphatidic acid has a glycerol backbone with two fatty acid chains 
    attached via ester bonds and one phosphate group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
        
    # Look for exactly 2 ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

    # Look for phosphate group (-O-P(=O)(O)(O))
    phosphate_pattern = Chem.MolFromSmarts("[OX2][PX4](=[OX1])([OX2H,OX1-])[OX2H,OX1-]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Check for fatty acid chains (long carbon chains attached to esters)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Missing fatty acid chains, found {len(fatty_acid_matches)}"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 8:
        return False, "Chains too short to be fatty acids"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    
    if c_count < 10:
        return False, "Too few carbons for phosphatidic acid"
    if o_count < 7:
        return False, "Must have at least 7 oxygens (2 esters + phosphate)"
    if p_count != 1:
        return False, "Must have exactly one phosphorus atom"

    # Check molecular weight - phosphatidic acids typically >400 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for phosphatidic acid"

    return True, "Contains glycerol backbone with 2 fatty acid chains and phosphate group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:16337',
                          'name': 'phosphatidic acid',
                          'definition': 'A derivative of glycerol in which one '
                                        'hydroxy group, commonly but not '
                                        'necessarily primary, is esterified '
                                        'with phosphoric acid and the other '
                                        'two are esterified with fatty acids.',
                          'parents': ['CHEBI:37739'],
                          'xrefs': ['KEGG:C00416'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'S(=O)(=O)(C1=C2C(=CC3=C1NC=4C=CC=CC34)[C@@]5([C@H]([C@](C(=O)O)([C@@H](O)CC5)C)CC2)C)C6=CC7=C(NC8=C7C=C9[C@@]%10([C@H]([C@](C(=O)O)([C@@H](O)CC%10)C)CCC9=C8)C)C=C6',
                                     'name': 'Sulfadixiamycin C',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'CNC(O)=O',
                                     'name': 'methylcarbamic acid',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'CCNC(=O)NC1=CC2=C(C=C1)OC[C@H]3[C@@H](CC[C@H](O3)CC(=O)N[C@@H](C)C4=CC=CC=C4)N(C2=O)C',
                                     'name': '2-[(2S,4aR,12aR)-8-(ethylcarbamoylamino)-5-methyl-6-oxo-2,3,4,4a,12,12a-hexahydropyrano[2,3-c][1,5]benzoxazocin-2-yl]-N-[(1S)-1-phenylethyl]acetamide',
                                     'reason': 'Found 0 ester groups, need '
                                               'exactly 2'},
                                 {   'smiles': 'O([C@H]1O[C@@H]([C@@H](O[C@@H]2O[C@@H]([C@@H](O[C@@H]3O[C@@H]([C@H](O)[C@H](O)[C@H]3O)CO[C@]4(O[C@H]([C@H](NC(=O)C)[C@@H](O)C4)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]2NC(=O)C)CO)[C@H](O)[C@@H]1O[C@@H]5O[C@@H]([C@@H](O[C@@H]6O[C@@H]([C@H](O)[C@H](O)[C@H]6O)CO[C@]7(O[C@H]([C@H](NC(=O)C)[C@@H](O)C7)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]5NC(=O)C)CO)CO)[C@H]8[C@H](O)[C@H](O[C@@H](O[C@H]9[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]9CO)O[C@H]%10[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]%10CO[C@@H]%11O[C@H]([C@@H](O)[C@@H](O)[C@@H]%11O)C)O)[C@H]8O)CO[C@H]%12O[C@@H]([C@@H](O)[C@H](O)[C@@H]%12O[C@@H]%13O[C@@H]([C@@H](O[C@@H]%14O[C@@H]([C@H](O)[C@H](O[C@]%15(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%15)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]%14O)CO)[C@H](O)[C@H]%13NC(=O)C)CO)CO[C@@H]%16O[C@@H]([C@@H](O[C@@H]%17O[C@@H]([C@H](O)[C@H](O)[C@H]%17O)CO[C@]%18(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%18)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]%16NC(=O)C)CO',
                                     'name': 'CID 91851985',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'O(C(=O)C(C1C(CN2C(C1)C=3NC=4C(C3CC2)=CC=CC4)CC)=COC)C',
                                     'name': 'Methyl '
                                             '2-(3-ethyl-1,2,3,4,6,7,12,12b-octahydroindolo[2,3-a]quinolizin-2-yl)-3-methoxyprop-2-enoate',
                                     'reason': 'Found 1 ester groups, need '
                                               'exactly 2'},
                                 {   'smiles': 'O[C@H](/C=C/C=C/C=C/[C@H](O)[C@H](O)C=C)[C@H](O)/C=C/C',
                                     'name': 'Separacene C',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'C[C@@H]1O[C@@H](O[C@@H]2[C@@H](CO)O[C@@H](O[C@@H]3[C@@H](O)C(O)O[C@H](CO)[C@@H]3O)[C@H](NC(C)=O)[C@H]2O[C@@H]2O[C@H](CO)[C@H](O)[C@H](OS(O)(=O)=O)[C@H]2O)[C@@H](O)[C@H](O)[C@@H]1O',
                                     'name': 'alpha-L-Fucp-(1->4)-[beta-D-Galp3S-(1->3)]-beta-D-GlcpNAc-(1->3)-D-Galp',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'C1=CC=CC2=C1C(N([C@H](C(N2)=O)CC=3C=CC(=CC3)OC)C)=O',
                                     'name': "(S)-4'-methoxycyclopeptine",
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'O=C(N[C@@H](CC=1C=2C(NC1)=CC=CC2)C(O)=O)[C@@H](NC(=O)[C@@H](N)C(C)C)C(C)C',
                                     'name': 'Val-Val-Trp',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'C=1C(=C(C=CC1/C=C/CO)OC(CO)C(O)C=2C=C(C(=CC2)O)OC)OC',
                                     'name': 'guaiacylglycerol beta-coniferyl '
                                             'ether',
                                     'reason': 'No glycerol backbone found'}],
    'sample_false_negatives': [   {   'smiles': 'CCCCC\\C=C/C\\C=C/C\\C=C/CCCCC(=O)OC[C@@H](O)COP(O)(O)=O',
                                      'name': '1-(gamma-linolenoyl)-sn-glycero-3-phosphate',
                                      'reason': 'Found 1 ester groups, need '
                                                'exactly 2'},
                                  {   'smiles': 'P(OC[C@H](OC(=O)CCCCCCCC(O)/C=C/C(O)=O)COC(=O)CCCCCCC/C=C\\CCCCCCCC)(O)(O)=O',
                                      'name': 'OHDdiA-PA',
                                      'reason': 'Found 3 ester groups, need '
                                                'exactly 2'},
                                  {   'smiles': 'CCCCCCCCCCCCCCCC(=O)OCC(COP(O)(=O)OCC(O)COP(O)(O)=O)OC(=O)CCCCCCC\\C=C/CCCCCCCC',
                                      'name': '1-hexadecanoyl-2-[(Z)-octadec-9-enoyl]-sn-glycero-3-phospho-sn-glycerol '
                                              '3-phosphate',
                                      'reason': 'Must have exactly one '
                                                'phosphorus atom'},
                                  {   'smiles': 'P(OC[C@@H](COC(CCCCC)=O)OC(CCCCC)=O)(=O)(O)O',
                                      'name': '1,2-dihexanoyl-sn-glycero-3-phosphate',
                                      'reason': 'Molecular weight too low for '
                                                'phosphatidic acid'},
                                  {   'smiles': 'P(OC[C@H](OC(=O)CCCCCCCC(=O)/C=C/C(O)=O)COC(=O)CCCCCCCCCCCCCCC)(O)(O)=O',
                                      'name': 'PKDdiA-PA',
                                      'reason': 'Found 3 ester groups, need '
                                                'exactly 2'},
                                  {   'smiles': 'O1C(C1)COC(=O)CCCC=CCC=CCC=CCC=CCCCCC',
                                      'name': 'O-Arachidonoyl Glycidol',
                                      'reason': 'Found 1 ester groups, need '
                                                'exactly 2'},
                                  {   'smiles': 'P(OC[C@H](OC(=O)CCC(O)/C=C/C(O)=O)COC(=O)CCCCCCC/C=C\\CCCCCCCC)(O)(O)=O',
                                      'name': 'OHHdiA-PA',
                                      'reason': 'Found 3 ester groups, need '
                                                'exactly 2'},
                                  {   'smiles': 'P(OC[C@H](OC(=O)CCC(=O)/C=C/C(O)=O)COC(=O)CCCCCCCCCCCCCCC)(O)(O)=O',
                                      'name': 'PKHdiA-PA',
                                      'reason': 'Found 3 ester groups, need '
                                                'exactly 2'},
                                  {   'smiles': 'P(OC[C@H](OC(=O)CCCCCCCC(O)/C=C/C(O)=O)COC(=O)CCCCCCCCCCCCCCC)(O)(O)=O',
                                      'name': 'PHDdiA-PA',
                                      'reason': 'Found 3 ester groups, need '
                                                'exactly 2'},
                                  {   'smiles': 'P(OC[C@H](OC(=O)CCCC(O)/C=C/C(O)=O)COC(=O)CCCCCCC/C=C\\CCCCCCCC)(O)(O)=O',
                                      'name': 'OHODiA-PA',
                                      'reason': 'Found 3 ester groups, need '
                                                'exactly 2'}],
    'attempt': 0,
    'reasoning': 'Let me break this down step by step:\n'
                 '\n'
                 '1. A phosphatidic acid has a glycerol backbone with:\n'
                 '   - Two fatty acid chains attached via ester bonds\n'
                 '   - One phosphate group attached to the remaining OH group\n'
                 '\n'
                 '2. Key structural features to check:\n'
                 '   - Glycerol backbone (C-C-C)\n'
                 '   - Two ester groups (-O-C(=O)-)\n'
                 '   - One phosphate group (-O-P(=O)(OH)2)\n'
                 '   - Long carbon chains (fatty acids)\n'
                 '\n'
                 '3. Common characteristics from examples:\n'
                 '   - Various fatty acid chain lengths\n'
                 '   - Can have saturated or unsaturated fatty acids\n'
                 '   - Specific stereochemistry at the glycerol carbon\n'
                 '\n'
                 "Here's the code:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 351,
    'num_false_positives': 78,
    'num_true_negatives': 141857,
    'num_false_negatives': 14,
    'num_negatives': None,
    'precision': 0.8181818181818182,
    'recall': 0.9616438356164384,
    'f1': 0.8841309823677582,
    'accuracy': 0.999353478566409,
    'negative_predictive_value': 0.9999013188037019}