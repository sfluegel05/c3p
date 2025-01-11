"""
Classifies: CHEBI:166828 saccharolipid
"""
"""
Classifies: saccharolipids - lipids containing carbohydrate moieties
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.
    Saccharolipids are lipids that contain carbohydrate moieties.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a saccharolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for minimum complexity - saccharolipids are large molecules
    if mol.GetNumAtoms() < 20:
        return False, "Molecule too small to be a saccharolipid"

    # Look for sugar rings (pyranose/furanose)
    sugar_pattern = Chem.MolFromSmarts("[CR1]1[CR1][CR1][CR1][CR1][OR1]1")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "No sugar rings found"

    # Look for long carbon chains (lipid part)
    lipid_pattern = Chem.MolFromSmarts("[CH2][CH2][CH2][CH2][CH2][CH2]")
    lipid_matches = mol.GetSubstructMatches(lipid_pattern)
    if not lipid_matches:
        return False, "No long carbon chains found"

    # Count key elements typical for saccharolipids
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    num_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if num_carbons < 12:
        return False, "Insufficient carbon atoms for saccharolipid"
    if num_oxygens < 6:
        return False, "Insufficient oxygen atoms for saccharolipid"

    # Look for ester/amide linkages (common in saccharolipids)
    ester_pattern = Chem.MolFromSmarts("[#6]-C(=O)-O-[#6]")
    amide_pattern = Chem.MolFromSmarts("[#6]-C(=O)-N-[#6]")
    
    has_ester = mol.HasSubstructMatch(ester_pattern)
    has_amide = mol.HasSubstructMatch(amide_pattern)
    
    if not (has_ester or has_amide):
        return False, "No ester or amide linkages found"

    # Check for hydroxyl groups (characteristic of sugars)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:
        return False, "Insufficient hydroxyl groups for sugar moiety"

    # Optional: Check for phosphate groups (common in many saccharolipids)
    phosphate_pattern = Chem.MolFromSmarts("[P](=O)([O-,OH])([O-,OH])")
    has_phosphate = mol.HasSubstructMatch(phosphate_pattern)

    # Calculate molecular weight - saccharolipids are typically large molecules
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for saccharolipid"

    # Construct reason string
    reason = f"Contains {len(sugar_matches)} sugar ring(s), lipid chains, "
    reason += f"{len(hydroxyl_matches)} hydroxyl groups, "
    reason += "and appropriate linkages between sugar and lipid components"
    if has_phosphate:
        reason += " with phosphate groups"

    return True, reason


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:166828',
                          'name': 'saccharolipid',
                          'definition': 'Lipids that contain a carbohydrate '
                                        'moiety.',
                          'parents': ['CHEBI:18059'],
                          'xrefs': ['Wikipedia:Saccharolipid'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'S(=O)(=O)(C1=C2C(=CC3=C1NC=4C=CC=CC34)[C@@]5([C@H]([C@](C(=O)O)([C@@H](O)CC5)C)CC2)C)C6=CC7=C(NC8=C7C=C9[C@@]%10([C@H]([C@](C(=O)O)([C@@H](O)CC%10)C)CCC9=C8)C)C=C6',
                                     'name': 'Sulfadixiamycin C',
                                     'reason': 'No sugar rings found'},
                                 {   'smiles': 'CNC(O)=O',
                                     'name': 'methylcarbamic acid',
                                     'reason': 'Molecule too small to be a '
                                               'saccharolipid'},
                                 {   'smiles': 'CCNC(=O)NC1=CC2=C(C=C1)OC[C@H]3[C@@H](CC[C@H](O3)CC(=O)N[C@@H](C)C4=CC=CC=C4)N(C2=O)C',
                                     'name': '2-[(2S,4aR,12aR)-8-(ethylcarbamoylamino)-5-methyl-6-oxo-2,3,4,4a,12,12a-hexahydropyrano[2,3-c][1,5]benzoxazocin-2-yl]-N-[(1S)-1-phenylethyl]acetamide',
                                     'reason': 'No sugar rings found'},
                                 {   'smiles': 'O([C@H]1O[C@@H]([C@@H](O[C@@H]2O[C@@H]([C@@H](O[C@@H]3O[C@@H]([C@H](O)[C@H](O)[C@H]3O)CO[C@]4(O[C@H]([C@H](NC(=O)C)[C@@H](O)C4)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]2NC(=O)C)CO)[C@H](O)[C@@H]1O[C@@H]5O[C@@H]([C@@H](O[C@@H]6O[C@@H]([C@H](O)[C@H](O)[C@H]6O)CO[C@]7(O[C@H]([C@H](NC(=O)C)[C@@H](O)C7)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]5NC(=O)C)CO)CO)[C@H]8[C@H](O)[C@H](O[C@@H](O[C@H]9[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]9CO)O[C@H]%10[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]%10CO[C@@H]%11O[C@H]([C@@H](O)[C@@H](O)[C@@H]%11O)C)O)[C@H]8O)CO[C@H]%12O[C@@H]([C@@H](O)[C@H](O)[C@@H]%12O[C@@H]%13O[C@@H]([C@@H](O[C@@H]%14O[C@@H]([C@H](O)[C@H](O[C@]%15(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%15)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]%14O)CO)[C@H](O)[C@H]%13NC(=O)C)CO)CO[C@@H]%16O[C@@H]([C@@H](O[C@@H]%17O[C@@H]([C@H](O)[C@H](O)[C@H]%17O)CO[C@]%18(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%18)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]%16NC(=O)C)CO',
                                     'name': 'CID 91851985',
                                     'reason': 'No long carbon chains found'},
                                 {   'smiles': 'O(C(=O)C(C1C(CN2C(C1)C=3NC=4C(C3CC2)=CC=CC4)CC)=COC)C',
                                     'name': 'Methyl '
                                             '2-(3-ethyl-1,2,3,4,6,7,12,12b-octahydroindolo[2,3-a]quinolizin-2-yl)-3-methoxyprop-2-enoate',
                                     'reason': 'No sugar rings found'},
                                 {   'smiles': 'O[C@H](/C=C/C=C/C=C/[C@H](O)[C@H](O)C=C)[C@H](O)/C=C/C',
                                     'name': 'Separacene C',
                                     'reason': 'Molecule too small to be a '
                                               'saccharolipid'},
                                 {   'smiles': 'C[C@@H]1O[C@@H](O[C@@H]2[C@@H](CO)O[C@@H](O[C@@H]3[C@@H](O)C(O)O[C@H](CO)[C@@H]3O)[C@H](NC(C)=O)[C@H]2O[C@@H]2O[C@H](CO)[C@H](O)[C@H](OS(O)(=O)=O)[C@H]2O)[C@@H](O)[C@H](O)[C@@H]1O',
                                     'name': 'alpha-L-Fucp-(1->4)-[beta-D-Galp3S-(1->3)]-beta-D-GlcpNAc-(1->3)-D-Galp',
                                     'reason': 'No long carbon chains found'},
                                 {   'smiles': 'C1=CC=CC2=C1C(N([C@H](C(N2)=O)CC=3C=CC(=CC3)OC)C)=O',
                                     'name': "(S)-4'-methoxycyclopeptine",
                                     'reason': 'No sugar rings found'},
                                 {   'smiles': 'O=C(N[C@@H](CC=1C=2C(NC1)=CC=CC2)C(O)=O)[C@@H](NC(=O)[C@@H](N)C(C)C)C(C)C',
                                     'name': 'Val-Val-Trp',
                                     'reason': 'No sugar rings found'},
                                 {   'smiles': 'C=1C(=C(C=CC1/C=C/CO)OC(CO)C(O)C=2C=C(C(=CC2)O)OC)OC',
                                     'name': 'guaiacylglycerol beta-coniferyl '
                                             'ether',
                                     'reason': 'No sugar rings found'}],
    'sample_false_negatives': [   {   'smiles': 'CC(=O)OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O',
                                      'name': '6-O-acetyl-beta-D-glucose',
                                      'reason': 'Molecule too small to be a '
                                                'saccharolipid'},
                                  {   'smiles': 'OC[C@@H](O)C(=C)C(=O)OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O',
                                      'name': '6-tuliposide B',
                                      'reason': 'No long carbon chains found'},
                                  {   'smiles': 'O1[C@@H]([C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O[C@@H]/2OC=C(C(\\C2=C/C)CC(OC[C@H]3O[C@@H](OCCC4=CC=C(O)C=C4)[C@H](O)[C@@H](O)[C@@H]3O)=O)C(OC)=O)CO',
                                      'name': 'Nuezhenide',
                                      'reason': 'No long carbon chains found'}],
    'attempt': 0,
    'reasoning': "I'll create a program to classify saccharolipids. Based on "
                 'the examples, saccharolipids are complex molecules that '
                 'contain both carbohydrate (sugar) moieties and lipid '
                 'components. Key characteristics include:\n'
                 '\n'
                 '1. Presence of sugar units (like glucose, galactose, KDO)\n'
                 '2. Presence of fatty acid chains\n'
                 '3. Often contain phosphate groups\n'
                 '4. Linkages between sugar and lipid components\n'
                 '\n'
                 "Here's the implementation:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 63,
    'num_false_positives': 100,
    'num_true_negatives': 12677,
    'num_false_negatives': 3,
    'num_negatives': None,
    'precision': 0.38650306748466257,
    'recall': 0.9545454545454546,
    'f1': 0.5502183406113537,
    'accuracy': 0.9919800669625477,
    'negative_predictive_value': 0.999763406940063}