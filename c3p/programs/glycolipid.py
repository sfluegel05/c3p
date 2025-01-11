"""
Classifies: CHEBI:33563 glycolipid
"""
"""
Classifies: CHEBI:24400 glycolipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid based on its SMILES string.
    A glycolipid is a molecule with a carbohydrate (glycosyl) part linked to a lipid part via a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycosidic bond pattern (C-O-C between carbohydrate and lipid)
    glycosidic_pattern = Chem.MolFromSmarts("[C,c][OX2][C,c]")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No glycosidic bond found"

    # Look for carbohydrate part (mono-, di-, or tri-saccharide)
    carbohydrate_pattern = Chem.MolFromSmarts("[C,c][OX2][C,c][OX2][C,c]")
    carbohydrate_matches = mol.GetSubstructMatches(carbohydrate_pattern)
    if len(carbohydrate_matches) == 0:
        return False, "No carbohydrate part found"

    # Look for lipid part (long carbon chain or sphingosine derivative)
    lipid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    lipid_matches = mol.GetSubstructMatches(lipid_pattern)
    if len(lipid_matches) == 0:
        return False, "No lipid part found"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Chains too short to be lipid"

    # Check molecular weight - glycolipids typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for glycolipid"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 20:
        return False, "Too few carbons for glycolipid"
    if o_count < 5:
        return False, "Too few oxygens for glycolipid"

    return True, "Contains carbohydrate part linked to lipid part via glycosidic bond"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33563',
                          'name': 'glycolipid',
                          'definition': 'Any member of class of '
                                        '1,2-di-O-acylglycerols joined at '
                                        'oxygen 3 by a glycosidic linkage to a '
                                        'carbohydrate part (usually a mono-, '
                                        'di- or tri-saccharide). Some '
                                        'substances classified as bacterial '
                                        'glycolipids have the sugar part '
                                        'acylated by one or more fatty acids '
                                        'and the glycerol part may be absent.',
                          'parents': ['CHEBI:35740'],
                          'xrefs': ['KEGG:C05005', 'Wikipedia:Glycolipids'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'OC[C@H]1O[C@H](C[C@@H]1O)N1C=NC2=C1N=CNC[C@H]2O',
                                     'name': 'pentostatin',
                                     'reason': 'No carbohydrate part found'},
                                 {   'smiles': 'O=C(O)[C@]1([C@H]2[C@@](OC=3C=C4C5=C(C=CC=C5)NC4=CC3CC2)(CC[C@@H]1O)C)C',
                                     'name': 'Oxiamycin',
                                     'reason': 'No carbohydrate part found'},
                                 {   'smiles': 'C1CCC(C1)CC#CC2=CC=C(C=C2)[C@@H]3[C@@H]4CN(CC(=O)N4[C@@H]3CO)C(=O)C5=CC=C(C=C5)F',
                                     'name': '(6R,7R,8S)-7-[4-(3-cyclopentylprop-1-ynyl)phenyl]-4-[(4-fluorophenyl)-oxomethyl]-8-(hydroxymethyl)-1,4-diazabicyclo[4.2.0]octan-2-one',
                                     'reason': 'No glycosidic bond found'},
                                 {   'smiles': 'S(=O)(C(SSCCC)CC)CCC',
                                     'name': 'Propyl 1-(propylsulfinyl)propyl '
                                             'disulfide',
                                     'reason': 'No glycosidic bond found'},
                                 {   'smiles': 'ClC=1C(=O)[C@@H]([C@@](O)(C/C=C\\CCCCC)C1)C[C@H](OC(=O)C)[C@@H](OC(=O)C)CCCC(OC)=O',
                                     'name': 'punaglandin 6',
                                     'reason': 'No carbohydrate part found'},
                                 {   'smiles': 'CCCCCCCCCCCCCCCCCCCCCCCC[C@H](O)C(=O)N[C@@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)[C@H](O)[C@@H](O)CCCCCCCCCCCCCC',
                                     'name': 'Ins-1-P-Cer(t18:0/2-OH-26:0)',
                                     'reason': 'No glycosidic bond found'},
                                 {   'smiles': 'CCC(C)(C)C(=O)OC1CC(C=C2C1[C@H]([C@H](C=C2)C)CC[C@@H]3CC(CC(=O)O3)O)C',
                                     'name': '2,2-dimethylbutanoic acid '
                                             '[(7S,8S)-8-[2-[(2R)-4-hydroxy-6-oxo-2-oxanyl]ethyl]-3,7-dimethyl-1,2,3,7,8,8a-hexahydronaphthalen-1-yl] '
                                             'ester',
                                     'reason': 'No carbohydrate part found'},
                                 {   'smiles': '[C@@H]1(C2=CC=CC=C2)[C@@H](N)C1',
                                     'name': '(1S,2R)-tranylcypromine',
                                     'reason': 'No glycosidic bond found'},
                                 {   'smiles': 'OC(=O)/C=C/C#CC#CC#N',
                                     'name': '(2E)-7-Cyanohept-2-en-4,6-diynoic '
                                             'acid',
                                     'reason': 'No glycosidic bond found'},
                                 {   'smiles': 'O=C(O)CC=1N=C2[C@H](NC(=O)C)[C@H](O)[C@@H]([C@H](N2C1)CO)O',
                                     'name': 'Nagstatin',
                                     'reason': 'No glycosidic bond found'}],
    'sample_false_negatives': [   {   'smiles': 'CC(C)=CCC\\C(C)=C\\CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/COP(O)(=O)O[C@@H]1O[C@H](CO)[C@@H](O)C1=O',
                                      'name': 'trans,octacis-decaprenylphospho-beta-D-erythro-pentofuranosid-2-ulose',
                                      'reason': 'No carbohydrate part found'},
                                  {   'smiles': 'CC(CCOP(O)(=O)O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(/C)CC\\C=C(/C)CCC=C(C)C',
                                      'name': 'beta-D-mannosyl '
                                              'C35-phosphodolichol',
                                      'reason': 'No carbohydrate part found'},
                                  {   'smiles': 'CCCCCCCCCCCCC\\C=C\\[C@@H](O)[C@@H](N)CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O',
                                      'name': 'psychosine',
                                      'reason': 'Molecular weight too low for '
                                                'glycolipid'},
                                  {   'smiles': 'CC(C)=CCC\\C(C)=C\\CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/COP(O)(=O)O[C@@H]1O[C@H](CO)[C@@H](O)[C@@H]1O',
                                      'name': 'trans,octacis-decaprenylphospho-beta-D-arabinofuranose',
                                      'reason': 'No carbohydrate part found'},
                                  {   'smiles': 'CC(C)=CCC\\C(C)=C\\CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/COP(O)(=O)O[C@@H]1O[C@H](COP(O)(O)=O)[C@@H](O)[C@H]1O',
                                      'name': 'trans,octacis-decaprenylphospho-beta-D-ribofuranose '
                                              '5-phosphate',
                                      'reason': 'No carbohydrate part found'},
                                  {   'smiles': 'P(O[C@@H]1OC([C@H](O)C(O)C1O)C(O)=O)(OC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC\\C=C(\\CC/C=C(/CCC=C(C)C)\\C)/C)/C)/C)/C)/C)/C)/C)/C)/C)/C)(O)=O',
                                      'name': 'Dodecaprenyl '
                                              'phosphate-galacturonic acid',
                                      'reason': 'No carbohydrate part found'},
                                  {   'smiles': 'CC(=O)N[C@@H]1[C@@H](O)[C@H](O)[C@@H](CO)O[C@@H]1OP(O)(=O)OP(O)(=O)OC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(/C)CC\\C=C(/C)CCC=C(C)C',
                                      'name': 'N-acetyl-alpha-D-glucosaminyl-1-diphospho-ditrans,polycis-undecaprenol',
                                      'reason': 'No carbohydrate part found'},
                                  {   'smiles': 'CC(=O)N[C@@H]1[C@@H](O)[C@H](O)[C@@H](CO)O[C@@H]1OP(O)(=O)OP(O)(=O)OC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(/C)CCC=C(C)C',
                                      'name': 'N-acetyl-alpha-D-glucosaminyl-1-diphospho-trans,polycis-decaprenol',
                                      'reason': 'No carbohydrate part found'},
                                  {   'smiles': 'O1[C@@H]([C@@H](O)[C@H](O)[C@@H](O)C1OCC(N)C(O)/C=C/CCCCCCCCCCCCC)CO',
                                      'name': '(3R,4S,5S,6R)-2-{[(4E)-2-amino-3-hydroxyoctadec-4-en-1-yl]oxy}-6-(hydroxymethyl)oxane-3,4,5-triol',
                                      'reason': 'Molecular weight too low for '
                                                'glycolipid'},
                                  {   'smiles': 'CC(=O)N[C@@H]1[C@@H](O)[C@@H](O)[C@@H](CO)O[C@@H]1OP(O)(=O)OP(O)(=O)OC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(/C)CC\\C=C(/C)CCC=C(C)C',
                                      'name': 'N-acetyl-alpha-D-galactosaminyl-1-diphospho-ditrans,polycis-undecaprenol',
                                      'reason': 'No carbohydrate part found'}],
    'attempt': 0,
    'reasoning': 'To classify glycolipids, we need to identify molecules that '
                 'contain a carbohydrate (glycosyl) part linked to a lipid '
                 'part. The carbohydrate part is typically a mono-, di-, or '
                 'tri-saccharide, and the lipid part can be a fatty acid or a '
                 'sphingosine derivative. The glycosidic linkage connects the '
                 'carbohydrate to the lipid.\n'
                 '\n'
                 'Here’s the Python program to classify glycolipids:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 554,
    'num_false_positives': 100,
    'num_true_negatives': 838,
    'num_false_negatives': 36,
    'num_negatives': None,
    'precision': 0.8470948012232415,
    'recall': 0.9389830508474576,
    'f1': 0.8906752411575563,
    'accuracy': 0.9109947643979057,
    'negative_predictive_value': 0.9588100686498856}