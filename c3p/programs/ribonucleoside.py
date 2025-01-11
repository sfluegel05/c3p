"""
Classifies: CHEBI:18254 ribonucleoside
"""
"""
Classifies: CHEBI:18254 ribonucleoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_ribonucleoside(smiles: str):
    """
    Determines if a molecule is a ribonucleoside based on its SMILES string.
    A ribonucleoside consists of a nucleobase attached to a D-ribose sugar.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ribonucleoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for ribose pattern (furan ring with multiple OH groups)
    # [#6] represents carbon, [#8] represents oxygen
    # The @ symbols ensure correct stereochemistry for D-ribose
    ribose_pattern = Chem.MolFromSmarts("[#6]1-[#8]-[#6]-[#6]-[#6]1")
    if not mol.HasSubstructMatch(ribose_pattern):
        return False, "No ribose sugar found"

    # Check for hydroxyl groups on ribose (typically 3)
    hydroxyl_pattern = Chem.MolFromSmarts("[#8H1]")
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if hydroxyl_matches < 2:  # Allow for some modification of OH groups
        return False, f"Insufficient hydroxyl groups ({hydroxyl_matches})"

    # Look for nucleobase patterns
    # Pyrimidine pattern (6-membered ring with nitrogens)
    pyrimidine_pattern = Chem.MolFromSmarts("c1[n]c[n]c1")
    # Purine pattern (fused 5,6 ring system with nitrogens)
    purine_pattern = Chem.MolFromSmarts("c1[n]c2[n]c[n]c2[n]1")
    
    has_pyrimidine = mol.HasSubstructMatch(pyrimidine_pattern)
    has_purine = mol.HasSubstructMatch(purine_pattern)
    
    if not (has_pyrimidine or has_purine):
        return False, "No nucleobase (pyrimidine or purine) found"

    # Check for N-glycosidic bond connecting nucleobase to ribose
    # This is a nitrogen connected to the anomeric carbon of ribose
    n_glycosidic_pattern = Chem.MolFromSmarts("[#7]-[#6]1-[#8]-[#6]-[#6]-[#6]1")
    if not mol.HasSubstructMatch(n_glycosidic_pattern):
        return False, "No N-glycosidic bond found"

    # Additional check for reasonable molecular size
    if mol.GetNumAtoms() < 15:  # Minimum size for a basic nucleoside
        return False, "Molecule too small to be a ribonucleoside"
    
    # Check if the molecule has a reasonable number of rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 2:  # Should have at least ribose + base
        return False, "Insufficient ring count"

    # Success case - molecule appears to be a ribonucleoside
    base_type = "purine" if has_purine else "pyrimidine"
    return True, f"Contains ribose sugar connected to {base_type} base via N-glycosidic bond"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:18254',
                          'name': 'ribonucleoside',
                          'definition': 'Any nucleoside where the sugar '
                                        'component is D-ribose.',
                          'parents': ['CHEBI:33838', 'CHEBI:47019'],
                          'xrefs': ['KEGG:C00911'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'C[C@@H](C1=CC=CC=C1)NC(=O)C[C@H]2CC[C@@H]([C@@H](O2)CO)NC(=O)CN3CCOCC3',
                                     'name': '2-[(2R,5S,6R)-6-(hydroxymethyl)-5-[[2-(4-morpholinyl)-1-oxoethyl]amino]-2-oxanyl]-N-[(1S)-1-phenylethyl]acetamide',
                                     'reason': 'No ribose sugar found'},
                                 {   'smiles': 'C[N+](C)(C)[C@@H](Cc1c[nH]c(n1)S(=O)C[C@H](NC(=O)CC[C@H]([NH3+])C([O-])=O)C([O-])=O)C([O-])=O',
                                     'name': 'N(alpha)-(L-gamma-glutamyl)-hercynyl-L-cysteine '
                                             'sulfoxide(1-)',
                                     'reason': 'No ribose sugar found'},
                                 {   'smiles': 'O(C1=CC=2[C@]3([C@](N(CC3)C)(N(C2C=C1)C)[H])C)C(=O)N4CCC=5C(C4)=CC=CC5',
                                     'name': 'quilostigmine',
                                     'reason': 'No ribose sugar found'},
                                 {   'smiles': 'O[C@@H]1[C@]23[C@@]4(N(C[C@@]([C@]2(C[C@@]4([C@]56[C@]3(CC(=O)[C@](C5)(C([C@H]6O)=C)[H])[H])[H])[H])(CC1)C)CC)[H]',
                                     'name': 'Bullatine G',
                                     'reason': 'No ribose sugar found'},
                                 {   'smiles': 'O1C2(C(C3(C(C4(C(CC3OC(=O)C)C(OC(=O)CC4)(C)C)C)CC2)C)CC15C6N(C=7C5=CC=CC7)C(=O)C(N6)C)C',
                                     'name': 'Teraspiridole C_130091',
                                     'reason': 'Insufficient hydroxyl groups '
                                               '(0)'},
                                 {   'smiles': 'O(C1[C@@H](OC(=O)C)C(O[C@@H](OC2=C(OC3=C(C2=O)C(O)=CC(O[C@@H]4OC([C@@H](O)[C@H](O)C4O)CO)=C3CC=C(C)C)C5=CC=C(OC)C=C5)[C@H]1O)C)[C@@H]6OC[C@@H](O)[C@H](OC(=O)C)C6O',
                                     'name': 'Sempervirenoside A',
                                     'reason': 'No ribose sugar found'},
                                 {   'smiles': 'COC[C@]1(C(=O)C2CCN1CC2)CO',
                                     'name': '(2S)-2-(hydroxymethyl)-2-(methoxymethyl)-1-azabicyclo[2.2.2]octan-3-one',
                                     'reason': 'No ribose sugar found'},
                                 {   'smiles': 'Oc1c(C2CC(Cc3ccccc23)c2ccc(OCc3ccc(cc3)C(F)(F)F)cc2)c(=O)oc2ccccc12',
                                     'name': 'Flocoumafen',
                                     'reason': 'No ribose sugar found'},
                                 {   'smiles': 'O[C@H]1CC=2C(N(C=3C1=CC=CC3)C(=O)N)=CC=CC2',
                                     'name': '(S)-MHD',
                                     'reason': 'No ribose sugar found'},
                                 {   'smiles': 'S(OC=1C(O)=C(\\C=C\\C2=CC(O)=C(CC=C(C)C)C(O)=C2)C=CC1O)(O)(=O)=O',
                                     'name': '3-{(e)-2-[3,5-dihydroxy-4-(3-methyl-2-buten-1-yl)phenyl]vinyl}-2,6-dihydroxyphenyl '
                                             'hydrogen sulfate',
                                     'reason': 'No ribose sugar found'}],
    'sample_false_negatives': [   {   'smiles': 'OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1c(cc(=O)[nH]c1=O)C(O)=O',
                                      'name': 'orotidine',
                                      'reason': 'No nucleobase (pyrimidine or '
                                                'purine) found'},
                                  {   'smiles': 'OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1ccc(=S)[nH]c1=O',
                                      'name': '4-thiouridine',
                                      'reason': 'No nucleobase (pyrimidine or '
                                                'purine) found'},
                                  {   'smiles': 'Nc1ccn([C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O)c(=O)n1',
                                      'name': 'cytidine',
                                      'reason': 'No nucleobase (pyrimidine or '
                                                'purine) found'},
                                  {   'smiles': 'COC(=O)Cc1cn([C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O)c(=S)[nH]c1=O',
                                      'name': '5-methoxycarbonylmethyl-2-thiouridine',
                                      'reason': 'No nucleobase (pyrimidine or '
                                                'purine) found'},
                                  {   'smiles': 'Nc1nc(=O)n(cc1F)[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O',
                                      'name': '5-fluorocytidine',
                                      'reason': 'No nucleobase (pyrimidine or '
                                                'purine) found'},
                                  {   'smiles': '[C@@H]1(N2C(NC(=O)C=C2)=O)O[C@H](COC(=O)C)[C@H]([C@H]1OC(=O)C)OC(=O)C',
                                      'name': 'uridine triacetate',
                                      'reason': 'Insufficient hydroxyl groups '
                                                '(0)'},
                                  {   'smiles': 'OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cc(CNCC(O)=O)c(=O)[nH]c1=S',
                                      'name': '5-carboxymethylaminomethyl-2-thiouridine',
                                      'reason': 'No nucleobase (pyrimidine or '
                                                'purine) found'},
                                  {   'smiles': 'OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)N1C=CCNC1=O',
                                      'name': '3,4-dihydrozebularine',
                                      'reason': 'No nucleobase (pyrimidine or '
                                                'purine) found'},
                                  {   'smiles': 'NC(CCN1C(=O)CCN([C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O)C1=O)C(O)=O',
                                      'name': '3-(3-amino-3-carboxypropyl)-5,6-dihydrouridine',
                                      'reason': 'No nucleobase (pyrimidine or '
                                                'purine) found'},
                                  {   'smiles': 'N1=CN=C2C(=C1N)C(=CN2[C@H]3[C@H](O)[C@H](O)[C@H](O3)CO)C#N',
                                      'name': 'toyocamycin',
                                      'reason': 'No nucleobase (pyrimidine or '
                                                'purine) found'}],
    'attempt': 0,
    'reasoning': "I'll create a program to classify ribonucleosides. A "
                 'ribonucleoside consists of a nucleobase (purine or '
                 'pyrimidine) attached to a D-ribose sugar. The key structural '
                 'features are:\n'
                 '\n'
                 '1. D-ribose sugar with 4 carbons and 1 oxygen in a furan '
                 'ring\n'
                 "2. Multiple hydroxyl groups on the ribose (typically at 2', "
                 "3', and 5' positions)\n"
                 '3. A nucleobase (purine or pyrimidine) connected via an '
                 'N-glycosidic bond\n'
                 '4. The connection point is typically N1 for pyrimidines and '
                 'N9 for purines\n'
                 '\n'
                 "Here's the code:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 79,
    'num_false_positives': 100,
    'num_true_negatives': 6268,
    'num_false_negatives': 45,
    'num_negatives': None,
    'precision': 0.441340782122905,
    'recall': 0.6370967741935484,
    'f1': 0.5214521452145214,
    'accuracy': 0.9776648182378311,
    'negative_predictive_value': 0.9928718517345161}