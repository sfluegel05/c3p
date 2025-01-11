"""
Classifies: CHEBI:71548 dihydroagarofuran sesquiterpenoid
"""
"""
Classifies: dihydroagarofuran sesquiterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
from rdkit.Chem import rdMolDescriptors

def is_dihydroagarofuran_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a dihydroagarofuran sesquiterpenoid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (True/False for classification, reason for the classification)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Core structure pattern for dihydroagarofuran skeleton
    # Represents the characteristic tricyclic system with specific connectivity
    core_pattern = Chem.MolFromSmarts("[C]1[C][C]2[C]3[C][C][C]([C])[C]3(O[C]1([C])[C])[C]2")
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Missing dihydroagarofuran core structure"

    # Count rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 3:
        return False, "Insufficient ring count for dihydroagarofuran skeleton"

    # Check for ester groups (typically present in these compounds)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = len(mol.GetSubstructMatches(ester_pattern))
    if ester_matches < 2:
        return False, "Insufficient ester groups"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Should have at least 15 carbons (sesquiterpenoid core)
    if c_count < 15:
        return False, "Insufficient carbon count for sesquiterpenoid"
    
    # Should have multiple oxygens due to ester groups
    if o_count < 4:
        return False, "Insufficient oxygen count"

    # Check for characteristic bridged ring system
    bridged_pattern = Chem.MolFromSmarts("[C]12[C][C][C]1[C][C]2")
    if not mol.HasSubstructMatch(bridged_pattern):
        return False, "Missing characteristic bridged ring system"

    # Check for typical molecular weight range
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 1000:
        return False, "Molecular weight outside typical range for dihydroagarofuran sesquiterpenoids"

    # Count sp3 carbons (should have several)
    sp3_c = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[CX4]")))
    if sp3_c < 8:
        return False, "Insufficient sp3 carbons for dihydroagarofuran skeleton"

    return True, "Contains dihydroagarofuran skeleton with characteristic features"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:71548',
                          'name': 'dihydroagarofuran sesquiterpenoid',
                          'definition': 'Any sesquiterpenoid with a '
                                        'dihydroagarofuran skeleton.',
                          'parents': ['CHEBI:26658'],
                          'xrefs': ['PMID:17898902'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'O([C@H](COC(=O)CCCCCCCCCCCCC)COC(=O)CCCCCCCCCCC/C=C\\C/C=C\\CCCCC)C(=O)CCCCCCCCCCCCC',
                                     'name': 'TG(14:0/14:0/22:2(13Z,16Z))',
                                     'reason': 'Missing dihydroagarofuran core '
                                               'structure'},
                                 {   'smiles': 'O([C@@H]1O[C@@H]([C@H](O)[C@H](O[C@H]2O[C@@H]([C@H](O)[C@H](O)[C@H]2NC(=O)C)CO)[C@H]1O[C@@H]3O[C@H]([C@@H](O)[C@@H](O)[C@@H]3O)C)CO)[C@H]4[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]4CO)OC[C@@H](O)[C@H](O)[C@H](O[C@@H]5O[C@@H]([C@H](O)[C@H](O)[C@H]5O)CO)[C@@H](NC(=O)C)CO',
                                     'name': 'N-[(2R,3R,4R,5R,6R)-2-[(2S,3R,4S,5S,6R)-2-[(2R,3S,4R,5R,6R)-5-Acetamido-6-[(2R,3S,4R,5S)-5-acetamido-2,3,6-trihydroxy-4-[(2R,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyhexoxy]-4-hydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-5-hydroxy-6-(hydroxymethyl)-3-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-4-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]acetamide',
                                     'reason': 'Missing dihydroagarofuran core '
                                               'structure'},
                                 {   'smiles': 'P12(O[C@@]3(NC(=O)C)[C@@](O1)(O)[C@H](O)[C@H](O[C@]3(O2)O)CO)=O',
                                     'name': 'N-Acetylglucosamine phosphate',
                                     'reason': 'Missing dihydroagarofuran core '
                                               'structure'},
                                 {   'smiles': 'OC1=C(N=CC=C1)C',
                                     'name': '3-hydroxyl-2-methylpyridine',
                                     'reason': 'Missing dihydroagarofuran core '
                                               'structure'},
                                 {   'smiles': '[H][C@@]1(C)O[C@@]1([H])C1=C[C@H](O)[C@@H](C)OC1=O',
                                     'name': 'aspyrone',
                                     'reason': 'Missing dihydroagarofuran core '
                                               'structure'},
                                 {   'smiles': 'O=C1C(O)=C(NC2=C1C=C(C)C(=C2)C)C(=O)O',
                                     'name': '3,4-dihydroxy-6,7-dimethyl-quinoline-2-carboxylic '
                                             'acid',
                                     'reason': 'Missing dihydroagarofuran core '
                                               'structure'},
                                 {   'smiles': 'O1[C@@]2(OC(C)(C)[C@H]([C@H]2O)C)[C@H]([C@@H]3[C@@]4([C@@]1(C5=C(C(=C(O)C=C5)C)C=C4)CC3)C)C',
                                     'name': 'Blazeispirol D',
                                     'reason': 'Missing dihydroagarofuran core '
                                               'structure'},
                                 {   'smiles': 'CCS(=O)(=O)N1CC2(C1)CN([C@@H](C3=C2C4=C(N3C)C=C(C=C4)OC)CO)C',
                                     'name': "[(1S)-1'-ethylsulfonyl-7-methoxy-2,9-dimethyl-1-spiro[1,3-dihydropyrido[3,4-b]indole-4,3'-azetidine]yl]methanol",
                                     'reason': 'Missing dihydroagarofuran core '
                                               'structure'},
                                 {   'smiles': 'CN(C)CC(=O)N1CCC2(CC1)CN([C@@H](C3=C2C4=C(N3)C=C(C=C4)OC)CO)C(=O)C5=CN=CC=C5',
                                     'name': "2-(dimethylamino)-1-[(1S)-1-(hydroxymethyl)-7-methoxy-2-[oxo(3-pyridinyl)methyl]-1'-spiro[3,9-dihydro-1H-pyrido[3,4-b]indole-4,4'-piperidine]yl]ethanone",
                                     'reason': 'Missing dihydroagarofuran core '
                                               'structure'},
                                 {   'smiles': 'O=C1C2(C(C(CC1(C(O)=C(C2=O)CC=C(C)C)CC=C(C)C)CC=C(C)C)(CCC=C(C)C)C)C(=O)C(CC)C',
                                     'name': 'Adhyperforin',
                                     'reason': 'Missing dihydroagarofuran core '
                                               'structure'}],
    'sample_false_negatives': [   {   'smiles': 'C1(C(O[C@@]2([C@](O)([C@@]34[C@H](OC(=O)C=5C=CC=CC5)[C@@]([H])([C@H]([C@H]([C@@]3([C@@H](OC(=O)C)[C@H]2OC(C)=O)COC(=O)C)OC(=O)C)OC(C)=O)[C@@](COC(C6=C(C1C)C=CN=C6)=O)(C)O4)C)[H])=O)C',
                                      'name': 'cangorinine E-1',
                                      'reason': 'Missing dihydroagarofuran '
                                                'core structure'},
                                  {   'smiles': '[H][C@]12C[C@H](OC(=O)c3ccccc3)[C@]3(C)[C@@H](OC(C)=O)[C@H](C[C@@H](C)[C@@]3(OC1(C)C)[C@@H]2OC(=O)c1ccoc1)OC(=O)c1ccoc1',
                                      'name': 'orbiculin F',
                                      'reason': 'Missing dihydroagarofuran '
                                                'core structure'},
                                  {   'smiles': '[C@@]1(C(O[C@@]2([C@](O)([C@@]34[C@H](OC(C)=O)[C@@]([H])([C@H]([C@H]([C@@]3([C@@H](OC(C=5C=CC=CC5)=O)[C@H]2OC(C)=O)COC(=O)C)OC(=O)C)OC(C)=O)[C@@](COC(C6=C([C@H]1C)C=CN=C6)=O)(C)O4)C)[H])=O)(C)O',
                                      'name': 'wilfordinine C',
                                      'reason': 'Missing dihydroagarofuran '
                                                'core structure'},
                                  {   'smiles': '[C@@H]1([C@H]([C@@]2([C@@H](OC(=O)C)[C@H]([C@]3([C@](O)([C@]42[C@H](OC(=O)C)[C@]1([H])[C@@](COC(C5=C(CC[C@@H](C(O3)=O)C)N=CC=C5)=O)(C)O4)C)[H])O)COC(=O)C)OC(=O)C)OC(C)=O',
                                      'name': '2-O-deacetyleuonine',
                                      'reason': 'Missing dihydroagarofuran '
                                                'core structure'},
                                  {   'smiles': 'C1(C(O[C@@]2([C@](O)([C@@]34[C@H](OC(=O)C=5C=CC=CC5)[C@@]([H])([C@H]([C@H]([C@@]3([C@@H](OC(=O)C)[C@H]2OC(C=6C=NC=CC6)=O)COC(=O)C)OC(=O)C)OC(C)=O)[C@@](COC(C7=C(C1C)N=CC=C7)=O)(C)O4)C)[H])=O)C',
                                      'name': 'hyponine D',
                                      'reason': 'Missing dihydroagarofuran '
                                                'core structure'},
                                  {   'smiles': '[C@@]12([C@@H]([C@@H]([C@@]3([C@H]([C@]14[C@]([C@H]([C@@H]([C@@H]2OC(=O)C)O)OC(C(CCC5=NC=CC=C5C(OC[C@@]3(O4)C)=O)(C)O)=O)(O)C)OC(C)=O)[H])OC(C)=O)OC(=O)C)COC(C)=O',
                                      'name': 'wilforidine',
                                      'reason': 'Missing dihydroagarofuran '
                                                'core structure'},
                                  {   'smiles': '[H][C@]12C[C@H](OC(=O)c3ccccc3)[C@]3(COC(C)=O)[C@@H](OC(C)=O)[C@H](C[C@@H](C)[C@@]3(OC1(C)C)[C@@H]2OC(=O)c1ccccc1)OC(=O)c1ccccc1',
                                      'name': '15-acetoxyorbiculin G',
                                      'reason': 'Missing dihydroagarofuran '
                                                'core structure'},
                                  {   'smiles': '[C@@H]1([C@H]([C@@]2([C@@H](OC(=O)C)[C@H]([C@]3([C@]([C@]42[C@H](OC(=O)C)[C@]1([H])[C@@](COC(=O)C5=C(N=CC=C5)CC[C@@H](C(O3)=O)C)(C)O4)(C)O)[H])OC(C)=O)COC(=O)C)OC(=O)C)OC(C)=O',
                                      'name': 'wilformine',
                                      'reason': 'Missing dihydroagarofuran '
                                                'core structure'},
                                  {   'smiles': '[C@@H]1(C(O[C@@]2([C@](O)([C@@]34[C@H](O)[C@@]([H])([C@H]([C@H]([C@@]3([C@@H](OC(=O)C)[C@H]2OC(=O)C)COC(=O)C)OC(=O)C)OC(=O)C)[C@@](COC(C5=C(C1C)C=CN=C5)=O)(C)O4)C)[H])=O)C',
                                      'name': 'tripfordine C',
                                      'reason': 'Missing dihydroagarofuran '
                                                'core structure'},
                                  {   'smiles': '[H]C1OC(=O)c2cccnc2[C@@]([H])(C)[C@@]([H])(C)C(=O)O[C@@]2([H])[C@H](OC(=O)c3ccc(=O)n(C)c3)[C@H](OC(=O)c3ccccc3)[C@@]3(COC(C)=O)[C@H](OC(C)=O)[C@@H](O)[C@]4([H])[C@@H](OC(C)=O)[C@]3(O[C@@]14C)[C@@]2(C)O',
                                      'name': 'emarginatine F',
                                      'reason': 'Missing dihydroagarofuran '
                                                'core structure'}],
    'attempt': 0,
    'reasoning': 'Let me help create a program to classify dihydroagarofuran '
                 'sesquiterpenoids. Based on the examples provided, these '
                 'compounds have some key characteristics:\n'
                 '\n'
                 '1. They share a common core structure (dihydroagarofuran '
                 'skeleton)\n'
                 '2. They typically have multiple ester groups (acetoxy, '
                 'benzoyloxy, etc.)\n'
                 '3. They have a specific carbon framework with characteristic '
                 'ring system\n'
                 '4. They usually contain 15 carbons in the core (as '
                 'sesquiterpenoids)\n'
                 '\n'
                 "Here's the program to classify these compounds:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 142274,
    'num_false_negatives': 26,
    'num_negatives': None,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': 0.9998172874209417,
    'negative_predictive_value': 0.9998172874209417}