"""
Classifies: CHEBI:17855 triglyceride
"""
"""
Classifies: CHEBI:17855 triglyceride
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_triglyceride(smiles: str):
    """
    Determines if a molecule is a triglyceride based on its SMILES string.
    A triglyceride has a glycerol backbone with three fatty acid chains attached via ester bonds.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (True if molecule is a triglyceride, reason for classification)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with 3 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
        
    # Look for 3 ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 3:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 3"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 20:
        return False, "Too few carbons for triglyceride"
    if o_count != 6:
        return False, "Must have exactly 6 oxygens (3 ester groups)"

    # Check for fatty acid chains (long carbon chains attached to esters)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 3:
        return False, f"Missing fatty acid chains, found {len(fatty_acid_matches)}"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Chains too short to be fatty acids"

    # Check molecular weight - triglycerides typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight ({mol_wt:.1f}) too low for triglyceride"

    return True, "Contains glycerol backbone with 3 fatty acid chains attached via ester bonds"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:17855',
                          'name': 'triglyceride',
                          'definition': 'Any glyceride resulting from the '
                                        'condensation of all three hydroxy '
                                        'groups of glycerol '
                                        '(propane-1,2,3-triol) with fatty '
                                        'acids.',
                          'parents': ['CHEBI:47778', 'CHEBI:76579'],
                          'xrefs': [   'KEGG:C00422',
                                       'LIPID_MAPS_class:LMGL0301',
                                       'PMID:2474544'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'C[C@@H](C1=CC=CC=C1)NC(=O)C[C@H]2CC[C@@H]([C@@H](O2)CO)NC(=O)CN3CCOCC3',
                                     'name': '2-[(2R,5S,6R)-6-(hydroxymethyl)-5-[[2-(4-morpholinyl)-1-oxoethyl]amino]-2-oxanyl]-N-[(1S)-1-phenylethyl]acetamide',
                                     'reason': 'Found 0 ester groups, need '
                                               'exactly 3'},
                                 {   'smiles': 'C[N+](C)(C)[C@@H](Cc1c[nH]c(n1)S(=O)C[C@H](NC(=O)CC[C@H]([NH3+])C([O-])=O)C([O-])=O)C([O-])=O',
                                     'name': 'N(alpha)-(L-gamma-glutamyl)-hercynyl-L-cysteine '
                                             'sulfoxide(1-)',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'O(C1=CC=2[C@]3([C@](N(CC3)C)(N(C2C=C1)C)[H])C)C(=O)N4CCC=5C(C4)=CC=CC5',
                                     'name': 'quilostigmine',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'O[C@@H]1[C@]23[C@@]4(N(C[C@@]([C@]2(C[C@@]4([C@]56[C@]3(CC(=O)[C@](C5)(C([C@H]6O)=C)[H])[H])[H])[H])(CC1)C)CC)[H]',
                                     'name': 'Bullatine G',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'NC1=NC=NC2=C1N=CN2[C@@H]3O[C@H](COP(=O)(O)O)[C@@H](OC(=O)[C@@H](N)CCC(O)=O)[C@H]3O',
                                     'name': "3'-L-glutamyl-AMP",
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'O1C2(C(C3(C(C4(C(CC3OC(=O)C)C(OC(=O)CC4)(C)C)C)CC2)C)CC15C6N(C=7C5=CC=CC7)C(=O)C(N6)C)C',
                                     'name': 'Teraspiridole C_130091',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'O(C1[C@@H](OC(=O)C)C(O[C@@H](OC2=C(OC3=C(C2=O)C(O)=CC(O[C@@H]4OC([C@@H](O)[C@H](O)C4O)CO)=C3CC=C(C)C)C5=CC=C(OC)C=C5)[C@H]1O)C)[C@@H]6OC[C@@H](O)[C@H](OC(=O)C)C6O',
                                     'name': 'Sempervirenoside A',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'COC[C@]1(C(=O)C2CCN1CC2)CO',
                                     'name': '(2S)-2-(hydroxymethyl)-2-(methoxymethyl)-1-azabicyclo[2.2.2]octan-3-one',
                                     'reason': 'Found 0 ester groups, need '
                                               'exactly 3'},
                                 {   'smiles': 'Oc1c(C2CC(Cc3ccccc23)c2ccc(OCc3ccc(cc3)C(F)(F)F)cc2)c(=O)oc2ccccc12',
                                     'name': 'Flocoumafen',
                                     'reason': 'Found 0 ester groups, need '
                                               'exactly 3'},
                                 {   'smiles': 'O[C@H]1CC=2C(N(C=3C1=CC=CC3)C(=O)N)=CC=CC2',
                                     'name': '(S)-MHD',
                                     'reason': 'No glycerol backbone found'}],
    'sample_false_negatives': [   {   'smiles': 'CCCCC\\C=C/C=C/C(=O)OCC=C=CCCCC(=O)OCC(COC(=O)CCCCCCC\\C=C/C\\C=C/C\\C=C/CC)OC(=O)CCCCCCC\\C=C/C\\C=C/C\\C=C/CC',
                                      'name': '1,2-Di-(9Z,12Z,15Z-octadecatrienoyl)-3-(8-(2E,4Z-decadienoyloxy)-5,6-octadienoyl)-sn-glycerol',
                                      'reason': 'Found 4 ester groups, need '
                                                'exactly 3'},
                                  {   'smiles': 'C(C(COC(=O)CCCCCCC/C=C\\C[C@@H](CCCCCC)O)OC(=O)CCCCCCC/C=C\\C[C@@H](CCCCCC)O)OC(=O)CCCCCCC/C=C\\C[C@@H](CCCCCC)O',
                                      'name': 'triricinolein',
                                      'reason': 'Must have exactly 6 oxygens '
                                                '(3 ester groups)'},
                                  {   'smiles': 'O(C(=O)CCCCC)C(COC(=O)CCCCCCCCCCCCC)COC(=O)CCCCC',
                                      'name': '1,2-dicaproyl-3-myristoylglycerol',
                                      'reason': 'Molecular weight (498.4) too '
                                                'low for triglyceride'},
                                  {   'smiles': 'CC(=O)OCC(COC(C)=O)OC(C)=O',
                                      'name': 'triacetin',
                                      'reason': 'Too few carbons for '
                                                'triglyceride'},
                                  {   'smiles': 'CCCC(=O)OCC(COC(=O)CCC)OC(=O)CCC',
                                      'name': 'tributyrin',
                                      'reason': 'Too few carbons for '
                                                'triglyceride'},
                                  {   'smiles': 'CCCCCC(=O)OCC(COC(=O)CCCCC)OC(=O)CCCCC',
                                      'name': 'tricaproin',
                                      'reason': 'Molecular weight (386.3) too '
                                                'low for triglyceride'},
                                  {   'smiles': 'CCCCCCCC(=O)OCC(COC(=O)CCCCCCC)OC(=O)CCCCCCC',
                                      'name': 'trioctanoin',
                                      'reason': 'Molecular weight (470.4) too '
                                                'low for triglyceride'},
                                  {   'smiles': 'O(C(=O)CCCCCCCCCCCCCCCCC)CC(OC(=O)C)COC(=O)C',
                                      'name': '2,3-Diacetoxypropyl stearate',
                                      'reason': 'Molecular weight (442.3) too '
                                                'low for triglyceride'},
                                  {   'smiles': 'O([C@H](COC(=O)CCCCCCCCC)COC(=O)CCCCCCC)C(=O)CCCCCCC',
                                      'name': 'TG(8:0/8:0/10:0)',
                                      'reason': 'Molecular weight (498.4) too '
                                                'low for triglyceride'},
                                  {   'smiles': 'O(C(COC(=O)/C=C/C1=CC(OC)=C(O)C(OC)=C1)COC(=O)/C=C/C2=CC(OC)=C(O)C(OC)=C2)C(=O)/C=C/C3=CC(OC)=C(O)C(OC)=C3',
                                      'name': 'Glycerol trisinapate',
                                      'reason': 'Must have exactly 6 oxygens '
                                                '(3 ester groups)'}],
    'attempt': 0,
    'reasoning': "I'll create a program to classify triglycerides based on "
                 'their SMILES string. Triglycerides have a glycerol backbone '
                 'with three fatty acid chains attached via ester bonds. Key '
                 'characteristics include:\n'
                 '\n'
                 '1. Glycerol backbone (C-C-C with 3 oxygens attached)\n'
                 '2. Three ester groups (-O-C(=O)-)\n'
                 '3. Three fatty acid chains (long carbon chains)\n'
                 '4. Molecular weight typically >500 Da\n'
                 '5. Exactly 6 oxygens (3 ester groups)\n'
                 '\n'
                 "Here's the code:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 1174,
    'num_false_positives': 2,
    'num_true_negatives': 141113,
    'num_false_negatives': 11,
    'num_negatives': None,
    'precision': 0.9982993197278912,
    'recall': 0.9907172995780591,
    'f1': 0.9944938585345193,
    'accuracy': 0.9999086437104708,
    'negative_predictive_value': 0.9999220543635384}