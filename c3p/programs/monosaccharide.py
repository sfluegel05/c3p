"""
Classifies: CHEBI:35381 monosaccharide
"""
"""
Classifies: CHEBI:35381 monosaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a monosaccharide based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a monosaccharide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbons - must have at least 3
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 3:
        return False, "Less than 3 carbon atoms"
        
    # Count oxygens - must have multiple oxygen atoms
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 2:
        return False, "Too few oxygen atoms"
    
    # Look for hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if hydroxyl_matches < 2:
        return False, "Insufficient hydroxyl groups"

    # Check for aldehyde group (for aldoses)
    aldehyde_pattern = Chem.MolFromSmarts("[CH,CH2]=O")
    has_aldehyde = mol.HasSubstructMatch(aldehyde_pattern)
    
    # Check for ketone group (for ketoses)
    ketone_pattern = Chem.MolFromSmarts("[#6]-C(=O)-[#6]")
    has_ketone = mol.HasSubstructMatch(ketone_pattern)
    
    # Check for hemiacetal/hemiketal group (for cyclic forms)
    hemiacetal_pattern = Chem.MolFromSmarts("[OH]C1[O,C][C,O][C,O][C,O][C,O]1")
    hemiketal_pattern = Chem.MolFromSmarts("[OH]C1[O,C][C,O][C,O][C,O]1")
    has_cyclic_form = mol.HasSubstructMatch(hemiacetal_pattern) or mol.HasSubstructMatch(hemiketal_pattern)
    
    # Must have either aldehyde/ketone (open form) or hemiacetal/hemiketal (cyclic form)
    if not (has_aldehyde or has_ketone or has_cyclic_form):
        return False, "No aldehyde, ketone, or cyclic sugar form found"
    
    # Check molecular weight - should typically be under 300 Da for single sugar unit
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 300 and not (carbon_count <= 7):  # Allow slightly larger molecules if they have typical sugar carbon count
        return False, "Molecular weight too high for single sugar unit"
    
    # Check for disaccharide linkage pattern
    disaccharide_pattern = Chem.MolFromSmarts("[OX2]([#6])[#6]-[OX2]-[#6]")
    if mol.HasSubstructMatch(disaccharide_pattern):
        # Further check if this is just part of a cyclic structure
        if len(mol.GetSubstructMatches(disaccharide_pattern)) > 1:
            return False, "Contains glycosidic bonds (appears to be oligosaccharide)"
    
    # Calculate ring count
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 1:
        return False, "Contains multiple rings - likely not a monosaccharide"
        
    # Additional check for carbon-hydroxyl ratio (most carbons should have an OH)
    if hydroxyl_matches < (carbon_count - 2) and not has_cyclic_form:
        return False, "Insufficient hydroxyl groups for carbon count"

    return True, "Contains required monosaccharide features (multiple hydroxyls and aldehyde/ketone/cyclic form)"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35381',
                          'name': 'monosaccharide',
                          'definition': 'Parent monosaccharides are '
                                        'polyhydroxy aldehydes '
                                        'H[CH(OH)]nC(=O)H or polyhydroxy '
                                        'ketones H-[CHOH]n-C(=O)[CHOH]m-H with '
                                        'three or more carbon atoms. The '
                                        "generic term 'monosaccharide' (as "
                                        'opposed to oligosaccharide or '
                                        'polysaccharide) denotes a single '
                                        'unit, without glycosidic connection '
                                        'to other such units. It includes '
                                        'aldoses, dialdoses, aldoketoses, '
                                        'ketoses and diketoses, as well as '
                                        'deoxy sugars, provided that the '
                                        'parent compound has a (potential) '
                                        'carbonyl group.',
                          'parents': ['CHEBI:16646'],
                          'xrefs': ['KEGG:C06698'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [],
    'sample_false_negatives': [   {   'smiles': 'O[C@@H]([C@H](O)[C@H](O)C)[C@@H](O)CO',
                                      'name': '6-Deoxyglucitol',
                                      'reason': 'No aldehyde, ketone, or '
                                                'cyclic sugar form found'},
                                  {   'smiles': '[H][C@@]1(OC(=O)C(O)=C1O)[C@@H](O)CO',
                                      'name': 'L-ascorbic acid',
                                      'reason': 'No aldehyde, ketone, or '
                                                'cyclic sugar form found'},
                                  {   'smiles': 'C1(=C([C@@H]([C@H](C(O1)O)O)O)O)C(O)=O',
                                      'name': '4,5-dehydro-D-glucuronic acid',
                                      'reason': 'No aldehyde, ketone, or '
                                                'cyclic sugar form found'},
                                  {   'smiles': 'O=C(N[C@H]1C(O)C(O)[C@@H](O[C@@H]([C@H](N)[C@H](O)CO)[C@@H](O)CN)OC1CO)CC',
                                      'name': 'Sorbistin A1',
                                      'reason': 'Molecular weight too high for '
                                                'single sugar unit'},
                                  {   'smiles': '[S-]C(=O)C(O)CCC',
                                      'name': '2-Hydroxy-4-methylthio-butyrate',
                                      'reason': 'Insufficient hydroxyl groups'},
                                  {   'smiles': 'O=C(N[C@H]1C(O)C(O)[C@@H](O[C@@H]([C@H](N)[C@H](O)CO)[C@@H](O)CN)OC1CO)C',
                                      'name': 'Sorbistin B',
                                      'reason': 'Molecular weight too high for '
                                                'single sugar unit'},
                                  {   'smiles': 'OC[C@H](CC(C(=O)[H])=O)O',
                                      'name': '(S)-4,5-dihydroxy-2-oxopentanal',
                                      'reason': 'Insufficient hydroxyl groups '
                                                'for carbon count'},
                                  {   'smiles': 'C[C@@H]1O[C@@H](O)[C@H](O[C@H]2OC(=C[C@H](O)[C@H]2O)C(O)=O)[C@H](O)[C@H]1O',
                                      'name': '2-O-(4-deoxy-beta-L-threo-hex-4-enopyranuronosyl)-alpha-L-rhamnopyranose',
                                      'reason': 'Molecular weight too high for '
                                                'single sugar unit'},
                                  {   'smiles': 'C(CCC(=O)NO)CC(=O)NN=C[C@H]([C@H]([C@@H]([C@@H](CO)O)O)O)O',
                                      'name': 'LSM-36764',
                                      'reason': 'No aldehyde, ketone, or '
                                                'cyclic sugar form found'},
                                  {   'smiles': 'C[C@H](O)[C@H](O)CC(=O)C(O)=O',
                                      'name': '2-dehydro-3-deoxy-L-rhamnonic '
                                              'acid',
                                      'reason': 'Insufficient hydroxyl groups '
                                                'for carbon count'}],
    'attempt': 0,
    'reasoning': "I'll create a program to classify monosaccharides. Let's "
                 'break this down into key characteristics:\n'
                 '\n'
                 '1. Must have 3 or more carbon atoms\n'
                 '2. Contains multiple hydroxyl groups\n'
                 '3. Contains either an aldehyde group (aldoses) or a ketone '
                 'group (ketoses)\n'
                 '4. Can be in open chain or cyclic form (pyranose/furanose)\n'
                 '5. Can include deoxy sugars (missing some OH groups)\n'
                 '6. Must be a single unit (no glycosidic bonds between '
                 'multiple sugar units)\n'
                 '\n'
                 "Here's the implementation:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 242,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 100,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 0.7076023391812866,
    'f1': 0.8287671232876712,
    'accuracy': 0.7076023391812866,
    'negative_predictive_value': 0.0}