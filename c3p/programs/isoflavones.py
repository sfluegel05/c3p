"""
Classifies: CHEBI:38757 isoflavones
"""
"""
Classifies: CHEBI:47919 isoflavones
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_isoflavones(smiles: str):
    """
    Determines if a molecule is an isoflavone based on its SMILES string.
    Isoflavones have a 3-aryl-1-benzopyran-4-one (3-phenylchromen-4-one) skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic checks
    if mol.GetNumAtoms() < 15:  # Minimum atoms for isoflavone core
        return False, "Too few atoms for isoflavone"

    # Core structure patterns
    
    # 1. Basic isoflavone core - more flexible pattern
    # Matches the basic chromone skeleton with position 3 for aryl attachment
    isoflavone_core = Chem.MolFromSmarts('[#6]1~[#6]~[#6](=O)~c2cccc[c,C]2O1')
    
    # 2. Alternative core pattern to catch variations
    isoflavone_core_alt = Chem.MolFromSmarts('O=C1C=COc2ccccc12')
    
    # 3. Pattern for 3-aryl connection
    aryl_connection = Chem.MolFromSmarts('[#6]1~[#6](~[#6](=O)~c2cccc[c,C]2O1)~[#6]3[#6]~[#6]~[#6]~[#6]~[#6]3')

    if not (mol.HasSubstructMatch(isoflavone_core) or mol.HasSubstructMatch(isoflavone_core_alt)):
        return False, "Missing basic isoflavone core structure"

    # Check ring count (should have at least 3 rings)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 3:
        return False, "Insufficient number of rings"

    # Verify presence of ketone group
    ketone_pattern = Chem.MolFromSmarts('C(=O)')
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "Missing ketone group"

    # Check for allowed atoms
    # Isoflavones typically contain C, H, O and some common substituents
    allowed_atoms = {6, 1, 8, 7, 9, 17, 35, 53}  # C,H,O,N,F,Cl,Br,I
    atom_nums = {atom.GetAtomicNum() for atom in mol.GetAtoms()}
    if not atom_nums.issubset(allowed_atoms):
        return False, "Contains disallowed atoms for isoflavone"

    # Additional checks for common features
    
    # Count oxygen atoms (should have at least 2 for the basic structure)
    o_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8])
    if o_count < 2:
        return False, "Insufficient oxygen atoms for isoflavone structure"

    # Check for aromatic rings (should have at least 2)
    aromatic_rings = 0
    for atom in mol.GetAtoms():
        if atom.GetIsAromatic():
            aromatic_rings += 1
    if aromatic_rings < 8:  # Minimum number of aromatic atoms in basic structure
        return False, "Insufficient aromatic character"

    # If we've passed all checks, this is likely an isoflavone
    return True, "Contains isoflavone core structure with appropriate substitution pattern"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:38757',
                          'name': 'isoflavones',
                          'definition': 'Any isoflavonoid with a '
                                        '3-aryl-1-benzopyran-4-one '
                                        '(3-aryl-4H-chromen-4-one) skeleton '
                                        'and its substituted derivatives.',
                          'parents': ['CHEBI:3992', 'CHEBI:50753'],
                          'xrefs': ['MetaCyc:Isoflavones'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: NONE\n'
               'False positives: NONE\n'
               'False negatives: SMILES: '
               'O1C(=C(C(=O)C2=C1C=C(O)C=C2)C3=CC=CC=C3)C NAME: '
               '7-Hydroxy-2-methylisoflavone REASON: MISSED Missing basic '
               'isoflavone core structure\n'
               ' * SMILES: O1C2=C(C(=O)C(C3=CC(=C(O)C=C3)C=O)=C1)C=CC(O)=C2 '
               'NAME: Corylinal REASON: MISSED Missing basic isoflavone core '
               'structure\n'
               ' * SMILES: Oc1ccc(cc1)-c1coc2cc(O)c(O)cc2c1=O NAME: '
               "4',6,7-trihydroxyisoflavone REASON: MISSED Missing basic "
               'isoflavone core structure\n'
               ' * SMILES: '
               'COc1c(N)c(O)c(CC=C(C)C)c(c1CC=C(C)C)-c1coc2cc(O)cc(O)c2c1=O '
               'NAME: piscerythramine REASON: MISSED Missing basic isoflavone '
               'core structure\n'
               ' * SMILES: Oc1ccc2c(c1)occ(-c1cc3OCOc3cc1O)c2=O NAME: '
               "2',7-dihydroxy-4',5'-methylenedioxyisoflavone REASON: MISSED "
               'Missing basic isoflavone core structure\n'
               ' * SMILES: COc1ccc(cc1)-c1coc2cc(O)c(OC)c(O)c2c1=O NAME: '
               'Irisolidone REASON: MISSED Missing basic isoflavone core '
               'structure\n'
               ' * SMILES: COC1=CC2=C(C=C1)C(=O)C(=CO2)C3=CC4=C(C=C3OC)OCO4 '
               'NAME: '
               '7-methoxy-3-(6-methoxy-1,3-benzodioxol-5-yl)-1-benzopyran-4-one '
               'REASON: MISSED Missing basic isoflavone core structure\n'
               ' * SMILES: CC(C)=CCc1c(O)ccc(c1O)-c1coc2cc(O)cc(O)c2c1=O NAME: '
               'licoisoflavone A REASON: MISSED Missing basic isoflavone core '
               'structure\n'
               ' * SMILES: '
               'O1C2=C(C(=O)C(C3=CC=C(OC)C=C3)=C1)C=CC(OCC=C(C)C)=C2 NAME: '
               'Maximaisoflavone J REASON: MISSED Missing basic isoflavone '
               'core structure\n'
               ' * SMILES: BrC1=C2OC=C(C3=CC(Br)=C(OC)C(=C3)Br)C(C2=CC=C1O)=O '
               'NAME: Maduraktermol J REASON: MISSED Missing basic isoflavone '
               'core structure\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [],
    'sample_false_negatives': [   {   'smiles': 'O1C(=C(C(=O)C2=C1C=C(O)C=C2)C3=CC=CC=C3)C',
                                      'name': '7-Hydroxy-2-methylisoflavone',
                                      'reason': 'Missing basic isoflavone core '
                                                'structure'},
                                  {   'smiles': 'O1C2=C(C(=O)C(C3=CC(=C(O)C=C3)C=O)=C1)C=CC(O)=C2',
                                      'name': 'Corylinal',
                                      'reason': 'Missing basic isoflavone core '
                                                'structure'},
                                  {   'smiles': 'Oc1ccc(cc1)-c1coc2cc(O)c(O)cc2c1=O',
                                      'name': "4',6,7-trihydroxyisoflavone",
                                      'reason': 'Missing basic isoflavone core '
                                                'structure'},
                                  {   'smiles': 'COc1c(N)c(O)c(CC=C(C)C)c(c1CC=C(C)C)-c1coc2cc(O)cc(O)c2c1=O',
                                      'name': 'piscerythramine',
                                      'reason': 'Missing basic isoflavone core '
                                                'structure'},
                                  {   'smiles': 'Oc1ccc2c(c1)occ(-c1cc3OCOc3cc1O)c2=O',
                                      'name': "2',7-dihydroxy-4',5'-methylenedioxyisoflavone",
                                      'reason': 'Missing basic isoflavone core '
                                                'structure'},
                                  {   'smiles': 'COc1ccc(cc1)-c1coc2cc(O)c(OC)c(O)c2c1=O',
                                      'name': 'Irisolidone',
                                      'reason': 'Missing basic isoflavone core '
                                                'structure'},
                                  {   'smiles': 'COC1=CC2=C(C=C1)C(=O)C(=CO2)C3=CC4=C(C=C3OC)OCO4',
                                      'name': '7-methoxy-3-(6-methoxy-1,3-benzodioxol-5-yl)-1-benzopyran-4-one',
                                      'reason': 'Missing basic isoflavone core '
                                                'structure'},
                                  {   'smiles': 'CC(C)=CCc1c(O)ccc(c1O)-c1coc2cc(O)cc(O)c2c1=O',
                                      'name': 'licoisoflavone A',
                                      'reason': 'Missing basic isoflavone core '
                                                'structure'},
                                  {   'smiles': 'O1C2=C(C(=O)C(C3=CC=C(OC)C=C3)=C1)C=CC(OCC=C(C)C)=C2',
                                      'name': 'Maximaisoflavone J',
                                      'reason': 'Missing basic isoflavone core '
                                                'structure'},
                                  {   'smiles': 'BrC1=C2OC=C(C3=CC(Br)=C(OC)C(=C3)Br)C(C2=CC=C1O)=O',
                                      'name': 'Maduraktermol J',
                                      'reason': 'Missing basic isoflavone core '
                                                'structure'}],
    'attempt': 3,
    'reasoning': 'Let me analyze the previous attempt and suggest '
                 'improvements:\n'
                 '\n'
                 '1. The core issue seems to be that the SMARTS patterns were '
                 "too restrictive and didn't match the variety of isoflavone "
                 'structures shown in the examples.\n'
                 '\n'
                 '2. Key observations from examples:\n'
                 '- They all share a basic 3-phenylchromen-4-one skeleton\n'
                 "- There's significant variation in substitution patterns "
                 '(hydroxy, methoxy, prenyl groups, etc.)\n'
                 '- Some have additional fused rings or complex substituents\n'
                 '- The aryl group at position 3 can have various substitution '
                 'patterns\n'
                 '\n'
                 '3. Improvements needed:\n'
                 '- Simplify the core SMARTS pattern to be more flexible\n'
                 '- Allow for both aromatic and non-aromatic forms\n'
                 '- Better handle variations in the aryl group substitution\n'
                 '- Account for common modifications like glycosylation\n'
                 '\n'
                 "Here's the improved code:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 1,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 100,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 0.009900990099009901,
    'f1': 0.0196078431372549,
    'accuracy': 0.009900990099009901,
    'negative_predictive_value': 0.0}