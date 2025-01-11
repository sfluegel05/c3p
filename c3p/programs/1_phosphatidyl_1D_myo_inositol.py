"""
Classifies: CHEBI:16749 1-phosphatidyl-1D-myo-inositol
"""
"""
Classifies: CHEBI:28874 1-phosphatidyl-1D-myo-inositol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_1_phosphatidyl_1D_myo_inositol(smiles: str):
    """
    Determines if a molecule is a 1-phosphatidyl-1D-myo-inositol based on its SMILES string.
    A 1-phosphatidyl-1D-myo-inositol consists of a glycerol backbone with two fatty acid chains,
    a phosphate group, and a 1D-myo-inositol moiety attached to the phosphate group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-phosphatidyl-1D-myo-inositol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with 2 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
        
    # Look for 2 ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

    # Check for fatty acid chains (long carbon chains attached to esters)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"

    # Look for phosphate group (-O-P(=O)(-O)-)
    phosphate_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX2])[OX2]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Look for 1D-myo-inositol moiety (6-membered ring with 6 hydroxyl groups)
    inositol_pattern = Chem.MolFromSmarts("[C@H]1([C@@H]([C@H]([C@@H]([C@H]([C@@H]1O)O)O)O)O)")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No 1D-myo-inositol moiety found"

    # Check that the phosphate group is attached to the inositol moiety
    phosphate_inositol_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX2])[OX2][C@H]1([C@@H]([C@H]([C@@H]([C@H]([C@@H]1O)O)O)O)O)")
    if not mol.HasSubstructMatch(phosphate_inositol_pattern):
        return False, "Phosphate group not attached to 1D-myo-inositol moiety"

    return True, "Contains glycerol backbone with two fatty acid chains, a phosphate group, and a 1D-myo-inositol moiety attached to the phosphate group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:16749',
                          'name': '1-phosphatidyl-1D-myo-inositol',
                          'definition': 'A phosphatidylinositol in which the '
                                        'inositol moiety is the 1D-myo isomer '
                                        'and the phosphatidyl group is located '
                                        'at its position 1.',
                          'parents': ['CHEBI:28874'],
                          'xrefs': ['KEGG:C01194', 'PMID:28600633'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'COC(=O)C1=CO[C@@H](O[C@@H]2O[C@H](CO[C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)[C@@H](O)[C@H](O)[C@H]2O)[C@H]2[C@@H]1CC=C2CO',
                                     'name': 'Genipin 1-beta-gentiobioside',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'O=[N+]([O-])C1=C2C(NC=C2C[C@@H]3N(C(=O)[C@](O)(CC4=CC(O)=C(O)C=C4)N(C3=O)C)C)=CC=C1',
                                     'name': '3,4-Dihydroxyphenyl-thaxtomin A',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'C[C@@H]1CN(C(=O)CCCN2C(=CN=N2)CO[C@@H]1CN(C)CC3=CC=CC=C3OC)[C@H](C)CO',
                                     'name': '(8R,9S)-6-[(2R)-1-hydroxypropan-2-yl]-9-[[(2-methoxyphenyl)methyl-methylamino]methyl]-8-methyl-10-oxa-1,6,14,15-tetrazabicyclo[10.3.0]pentadeca-12,14-dien-5-one',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'COc1ccc(C[C@@H](C)[C@@H](C)Cc2ccc(O)c(OC)c2)cc1OC',
                                     'name': "(-)-(8R,8'S)-3,3',4-trimethoxy-4'-hydroxylignan",
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'O1C=C([C@H](O)[C@H]([C@H]1O)O)C(O)C=2OC=3C(O)=C(O)C=CC3C2',
                                     'name': 'Deuteromycol A',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': '[H]C1C[C@H](C)OC(=O)C[C@@H](NC(=O)[C@@H](Cc2c(Br)[nH]c3ccccc23)N(C)C(=O)[C@H](C)NC(=O)[C@@H](C)C\\C(C)=C\\1)c1ccc(O)cc1',
                                     'name': 'jaspamide H',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'C[C@@H]1CN(C(=O)C2=C(C(=CC=C2)NC(=O)C3CCCCC3)O[C@H]1CN(C)CC4=CC=C(C=C4)OC)[C@H](C)CO',
                                     'name': 'N-[(2R,3R)-5-[(2R)-1-hydroxypropan-2-yl]-2-[[(4-methoxyphenyl)methyl-methylamino]methyl]-3-methyl-6-oxo-3,4-dihydro-2H-1,5-benzoxazocin-10-yl]cyclohexanecarboxamide',
                                     'reason': 'Found 0 ester groups, need '
                                               'exactly 2'},
                                 {   'smiles': 'O[C@H](CCCCCCCCCCCCC)COC[C@@H](O)CO',
                                     'name': '1-O-(2R-hydroxy-pentadecyl)-sn-glycerol',
                                     'reason': 'Found 0 ester groups, need '
                                               'exactly 2'},
                                 {   'smiles': 'C[C@@H]1CN(C(=O)CC2=C(C=CC(=C2)NC(=O)NC3=CC=CC4=CC=CC=C43)O[C@@H]1CN(C)C(=O)NC5=CC=C(C=C5)C(F)(F)F)[C@@H](C)CO',
                                     'name': '1-[[(2S,3R)-5-[(2S)-1-hydroxypropan-2-yl]-3-methyl-9-[[(1-naphthalenylamino)-oxomethyl]amino]-6-oxo-2,3,4,7-tetrahydro-1,5-benzoxazonin-2-yl]methyl]-1-methyl-3-[4-(trifluoromethyl)phenyl]urea',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'O1[C@H]([C@H](NC(=O)C)[C@@H](O)C[C@]1(OC[C@H]2O[C@@H](O)[C@H](NC(=O)C)[C@@H](O[C@@H]3O[C@@H]([C@H](O)[C@H](O)[C@H]3O)CO)[C@@H]2O)C(O)=O)[C@H](O)[C@H](O)CO',
                                     'name': '(2R,4S,5R,6R)-5-Acetamido-2-[[(2R,3S,4R,5R,6R)-5-acetamido-3,6-dihydroxy-4-[(2R,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]methoxy]-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'No glycerol backbone found'}],
    'sample_false_negatives': [],
    'attempt': 0,
    'reasoning': 'To classify chemical entities as '
                 '1-phosphatidyl-1D-myo-inositol, we need to identify the key '
                 'structural features of this class. The molecule should '
                 'contain:\n'
                 '1. A glycerol backbone.\n'
                 '2. Two fatty acid chains attached to the glycerol via ester '
                 'bonds.\n'
                 '3. A phosphate group attached to the glycerol.\n'
                 '4. A 1D-myo-inositol moiety attached to the phosphate '
                 'group.\n'
                 '\n'
                 'The 1D-myo-inositol moiety is a specific stereoisomer of '
                 'inositol, which has six hydroxyl groups in a specific '
                 'spatial arrangement. The phosphate group is attached to the '
                 '1-position of the inositol ring.\n'
                 '\n'
                 'Here is the Python program to classify molecules as '
                 '1-phosphatidyl-1D-myo-inositol:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 27,
    'num_false_positives': 100,
    'num_true_negatives': 49901,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.2125984251968504,
    'recall': 1.0,
    'f1': 0.35064935064935066,
    'accuracy': 0.9980011193731511,
    'negative_predictive_value': 1.0}