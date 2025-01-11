"""
Classifies: CHEBI:50699 oligosaccharide
"""
"""
Classifies: oligosaccharide
A compound in which monosaccharide units are joined by glycosidic linkages.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_oligosaccharide(smiles: str):
    """
    Determines if a molecule is an oligosaccharide based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an oligosaccharide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for sugar rings with proper substitution patterns
    # Pyranose (6-membered) ring pattern with typical sugar substitutions
    pyranose_pattern = Chem.MolFromSmarts("[C]1[C]([OH0,OH1])[C]([OH0,OH1])[C]([OH0,OH1])[C]([OH0,OH1])O1")
    # Furanose (5-membered) ring pattern with typical sugar substitutions
    furanose_pattern = Chem.MolFromSmarts("[C]1[C]([OH0,OH1])[C]([OH0,OH1])[C]([OH0,OH1])O1")
    
    pyranose_matches = len(mol.GetSubstructMatches(pyranose_pattern))
    furanose_matches = len(mol.GetSubstructMatches(furanose_pattern))
    total_rings = pyranose_matches + furanose_matches
    
    if total_rings < 2:
        return False, f"Found only {total_rings} sugar rings, need at least 2"
    
    if total_rings > 20:
        return False, f"Found {total_rings} sugar rings, likely a polysaccharide"
        
    # Look for glycosidic linkages between rings
    # More specific pattern that looks for C-O-C between ring carbons
    glycosidic_pattern = Chem.MolFromSmarts("[C;R][O][C;R]")
    glycosidic_matches = len(mol.GetSubstructMatches(glycosidic_pattern))
    
    if glycosidic_matches < 1:
        return False, "No glycosidic linkages found"
        
    # Check for hydroxyl and hydroxyl-derived groups (including modified ones)
    hydroxyl_pattern = Chem.MolFromSmarts("[$([OH1]),$([OH0][C]=[O]),$(O[CH3]),$(O[C]=[O])]")
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
    
    if hydroxyl_matches < 4:
        return False, "Too few hydroxyl or hydroxyl-derived groups"
        
    # Count oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 6:
        return False, "Too few oxygen atoms for an oligosaccharide"
        
    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:  # Lowered threshold
        return False, "Molecular weight too low for oligosaccharide"
    if mol_wt > 10000:  # Increased threshold
        return False, "Molecular weight too high, likely a polysaccharide"
        
    # Check for characteristic sugar features
    # Look for CH2OH groups (common in sugars)
    ch2oh_pattern = Chem.MolFromSmarts("[CH2][OH1]")
    # Look for CHOH groups
    choh_pattern = Chem.MolFromSmarts("[CH1][OH1]")
    # Look for modified versions (acetylated, methylated, etc)
    modified_pattern = Chem.MolFromSmarts("[CH1,CH2][O][C]")
    
    sugar_features = (
        len(mol.GetSubstructMatches(ch2oh_pattern)) +
        len(mol.GetSubstructMatches(choh_pattern)) +
        len(mol.GetSubstructMatches(modified_pattern))
    )
    
    if sugar_features < 3:
        return False, "Missing characteristic sugar structural features"

    # Check for reasonable C:O ratio for sugars
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count > 0 and o_count/c_count < 0.3:
        return False, "Carbon to oxygen ratio not consistent with oligosaccharides"

    return True, f"Contains {total_rings} sugar rings connected by glycosidic linkages with appropriate sugar characteristics"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:50699',
                          'name': 'oligosaccharide',
                          'definition': 'A compound in which monosaccharide '
                                        'units are joined by glycosidic '
                                        'linkages. The term is commonly used '
                                        'to refer to a defined structure as '
                                        'opposed to a polymer of unspecified '
                                        'length or a homologous mixture. When '
                                        'the linkages are of other types the '
                                        'compounds are regarded as '
                                        'oligosaccharide analogues.',
                          'parents': ['CHEBI:16646', 'CHEBI:167559'],
                          'xrefs': ['KEGG:C00930'],
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
               'O1[C@@H]([C@H](O)[C@H](O)[C@@H](O)[C@H]1OC[C@H]2O[C@H](OC[C@H]3OC(O)[C@H](O)[C@@H](O)[C@@H]3O)[C@H](O)[C@@H](O)[C@H]2O)CO[C@@H]4O[C@@H]([C@H](O)[C@H](O)[C@H]4O)CO '
               'NAME: '
               '(3R,4S,5S,6R)-6-[[(2S,3R,4S,5R,6R)-3,4,5-Trihydroxy-6-[[(2S,3R,4S,5R,6R)-3,4,5-trihydroxy-6-[[(2R,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]oxan-2-yl]oxymethyl]oxan-2-yl]oxymethyl]oxane-2,3,4,5-tetrol '
               'REASON: MISSED Too few hydroxyl groups for an oligosaccharide\n'
               ' * SMILES: '
               'O1[C@@H]([C@@H](O)[C@H](O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)[C@@H](O)[C@@H]1OC[C@H]3O[C@@H](OC[C@H]4OC(O)[C@H](O)[C@@H](O[C@@H]5O[C@@H]([C@@H](O)[C@H](O)[C@H]5O)CO)[C@@H]4O)[C@H](O)[C@@H](O)[C@@H]3O)CO[C@@H]6O[C@@H]([C@@H](O)[C@H](O)[C@H]6O)CO '
               'NAME: '
               '(2R,3R,4S,5S,6R)-2-[[(2R,3R,4S,5R,6R)-3,5-Dihydroxy-4-[(2S,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-6-[[(2R,3S,4S,5R,6R)-3,4,5-trihydroxy-6-[[(2R,3R,4S,5R)-3,5,6-trihydroxy-4-[(2S,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]methoxy]oxan-2-yl]methoxy]oxan-2-yl]methoxy]-6-(hydroxymethyl)oxane-3,4,5-triol '
               'REASON: MISSED Too few hydroxyl groups for an oligosaccharide\n'
               ' * SMILES: '
               'O([C@@H]1[C@H](O)[C@@H](O)[C@@H](O[C@@H]1CO)O[C@H]2[C@H](O)[C@@H](O)[C@@H](OC2)O)[C@@H]3O[C@@H]([C@H](O)[C@H](O)[C@H]3O)CO '
               'NAME: '
               '(2S,3R,4S,5R,6R)-2-[(2R,3R,4R,5R,6R)-4,5-Dihydroxy-2-(hydroxymethyl)-6-[(3R,4R,5R,6R)-4,5,6-trihydroxyoxan-3-yl]oxyoxan-3-yl]oxy-6-(hydroxymethyl)oxane-3,4,5-triol '
               'REASON: MISSED Too few hydroxyl groups for an oligosaccharide\n'
               ' * SMILES: '
               'O([C@H]1O[C@@H]([C@@H](O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2NC(=O)C)CO)[C@H](O)[C@@H]1O[C@@H]3O[C@@H]([C@@H](O)[C@H](O)[C@H]3NC(=O)C)CO)CO)[C@H]4[C@H](O)[C@H](O[C@@H](O[C@@H]([C@H](O)[C@@H](NC(=O)C)CO)[C@H](O)CO)[C@H]4O)CO[C@H]5O[C@@H]([C@@H](O)[C@H](O)[C@@H]5O[C@@H]6O[C@@H]([C@@H](O)[C@H](O)[C@H]6NC(=O)C)CO)CO '
               'NAME: '
               'N-[(2S,3R,4R,5S,6R)-2-[(2R,3S,4S,5S,6R)-5-[(2S,3R,4R,5S,6R)-3-Acetamido-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-6-[(2R,3R,4S,5S,6S)-2-[[(2S,3S,4S,5S,6R)-3-[(2S,3R,4R,5S,6R)-3-acetamido-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]-6-[(2R,3S,4R,5S)-5-acetamido-1,2,4,6-tetrahydroxyhexan-3-yl]oxy-3,5-dihydroxyoxan-4-yl]oxy-4-hydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]acetamide '
               'REASON: MISSED Too few hydroxyl groups for an oligosaccharide\n'
               ' * SMILES: '
               'O=C(OC(CC1=CC=C(O[C@H]2OC([C@@H](OC)C(C2OC)O[C@H]3OC([C@@H](O)C(C3OC)O[C@H]4OC([C@H](OC(=O)C)C(C4OC)O[C@@H]5OC([C@@H](O)C(C5OC)O)CO)C)C)C)C=C1)CC(OC(=O)C(CCC)C)CCCCC(C(OC)CC)C)C(CCC)C '
               'NAME: PGL KIV REASON: MISSED Too few hydroxyl groups for an '
               'oligosaccharide\n'
               ' * SMILES: '
               'O([C@@H]1[C@@H](O)[C@H](O)[C@H](O[C@@H]1O)CO[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)[C@@H]3O[C@@H]([C@@H](O)[C@H](O)[C@H]3O)CO '
               'NAME: '
               '(2R,3S,4S,5R,6R)-2-(Hydroxymethyl)-6-[[(2R,3S,4S,5R,6S)-3,4,6-trihydroxy-5-[(2S,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]methoxy]oxane-3,4,5-triol '
               'REASON: MISSED Too few hydroxyl groups for an oligosaccharide\n'
               ' * SMILES: '
               'O1[C@@H]([C@H](O)[C@H](O)[C@@H](O)[C@H]1OC[C@H]2OC(O)[C@H](O)[C@@H](O)[C@H]2O)CO[C@H]3O[C@@H]([C@H](O)[C@H](O)[C@H]3O)CO '
               'NAME: '
               '(3R,4S,5R,6R)-6-[[(2S,3R,4S,5R,6R)-3,4,5-Trihydroxy-6-[[(2S,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]oxan-2-yl]oxymethyl]oxane-2,3,4,5-tetrol '
               'REASON: MISSED Too few hydroxyl groups for an oligosaccharide\n'
               ' * SMILES: '
               'O([C@H]1[C@H](O)[C@H](O[C@@H](O[C@H]2[C@H](O)[C@@H](O)C(O[C@@H]2CO)O)[C@H]1O)CO[C@H]3O[C@@H]([C@@H](O)[C@H](O)[C@@H]3O[C@@H]4O[C@@H]([C@@H](O[C@@H]5O[C@@H]([C@H](O)[C@H](O)[C@H]5O)CO)[C@H](O)[C@H]4NC(=O)C)CO)CO[C@@H]6O[C@@H]([C@@H](O[C@@H]7O[C@@H]([C@H](O)[C@H](O)[C@H]7O)CO)[C@H](O)[C@H]6NC(=O)C)CO)[C@H]8O[C@@H]([C@@H](O)[C@H](O)[C@@H]8O[C@@H]9O[C@@H]([C@@H](O[C@@H]%10O[C@@H]([C@H](O)[C@H](O)[C@H]%10O)CO)[C@H](O)[C@H]9NC(=O)C)CO)CO '
               'NAME: '
               'N-[(2R,3R,4R,5S,6R)-2-[[(2R,3S,4S,5S,6S)-5-[(2S,3R,4R,5S,6R)-3-Acetamido-4-hydroxy-6-(hydroxymethyl)-5-[(2S,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxy-6-[[(2R,3R,4S,5S,6S)-4-[(2R,3S,4S,5S,6R)-3-[(2S,3R,4R,5S,6R)-3-acetamido-4-hydroxy-6-(hydroxymethyl)-5-[(2S,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-3,5-dihydroxy-6-[(2R,3S,4R,5R)-4,5,6-trihydroxy-2-(hydroxymethyl)oxan-3-yl]oxyoxan-2-yl]methoxy]-3,4-dihydroxyoxan-2-yl]methoxy]-4-hydroxy-6-(hydroxymethyl)-5-[(2S,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-3-yl]acetamide '
               'REASON: MISSED Too few hydroxyl groups for an oligosaccharide\n'
               ' * SMILES: '
               '[H][C@@]1(O[C@@H]2O[C@H](COC(C)=O)[C@@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@H](O)[C@@H](COC(C)=O)O[C@H]1O[C@@]1([H])[C@@H](O[C@H](COC(C)=O)[C@@H](OC(C)=O)[C@@H]1OC(C)=O)O[C@]1([H])[C@@H](COC(C)=O)O[C@@H](O[C@]2([H])[C@@H](CO)O[C@@H](O[C@H](CCC)CCCCCCCCCCCCCCCCOC[C@@H](OC(C)=O)C(O)=O)[C@]([H])(O[C@@H]3O[C@H](COC(C)=O)[C@@H](OC(C)=O)[C@H](OC(C)=O)[C@H]3OC(C)=O)[C@H]2O)[C@H](OC(C)=O)[C@H]1OC(C)=O '
               'NAME: Agminoside A REASON: MISSED Too few hydroxyl groups for '
               'an oligosaccharide\n'
               ' * SMILES: '
               'O=C1O[C@@H]([C@@H](O[C@H]2O[C@@H]([C@@H](O)[C@@H]([C@H]2O)O)CO)[C@@H]([C@H]1O)O)CO '
               'NAME: D-maltobiono-1,5-lactone REASON: MISSED Too few hydroxyl '
               'groups for an oligosaccharide\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [],
    'sample_false_negatives': [   {   'smiles': 'O=C1O[C@@H]([C@@H](O[C@H]2O[C@@H]([C@@H](O)[C@@H]([C@H]2O)O)CO)[C@@H]([C@H]1O)O)CO',
                                      'name': 'D-maltobiono-1,5-lactone',
                                      'reason': 'Found only 1 sugar rings, '
                                                'need at least 2'},
                                  {   'smiles': '[H][C@@]1(O[C@@](C[C@@H](O)[C@H]1O)(OC[C@@H](O)[C@@]1([H])O[C@@](C[C@@H](O)[C@H]1O)(O[C@@H]1C[C@@](O)(O[C@]([H])([C@H](O)CO)[C@@H]1O)C(O)=O)C(O)=O)C(O)=O)[C@H](O)CO',
                                      'name': 'alpha-Kdo-(2->8)-alpha-Kdo-(2->4)-alpha-Kdo',
                                      'reason': 'Found only 0 sugar rings, '
                                                'need at least 2'},
                                  {   'smiles': 'O([C@H]1[C@H](O)[C@H](O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2NC(=O)C)CO)C(O[C@@H]1CO)O)[C@@H]3O[C@@H]([C@@H](O)[C@H](O)[C@H]3NC(=O)C)CO',
                                      'name': 'N-[(2S,3R,4R,5S,6R)-2-[(2R,3S,4S,5S)-5-[(2S,3R,4R,5S,6R)-3-Acetamido-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-4,6-dihydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]acetamide',
                                      'reason': 'Found only 1 sugar rings, '
                                                'need at least 2'},
                                  {   'smiles': 'O=C1C2=C(O)C3=C(O)C(=C(O[C@@H]4O[C@@H]([C@@H](O)[C@@H](C4)O)C)C=C3C=C2CC(C1O[C@@H]5O[C@@H]([C@@H](O)[C@@H](C5)O[C@@H]6O[C@@H]([C@H](O)[C@@H](C6)O[C@@H]7O[C@@H]([C@@H](OC)[C@@H](C7)O[C@@H]8O[C@@H]([C@@H](O)[C@@H](C8)O)C)C)C)C)C(OC)C(=O)C(O)C(O)C)C(CC)C',
                                      'name': 'UCH9',
                                      'reason': 'Found only 0 sugar rings, '
                                                'need at least 2'},
                                  {   'smiles': 'O([C@H]1[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]1CO)OC[C@@H](O)[C@H](O)[C@H](O)[C@@H](NC(=O)C)CO)[C@@H]2O[C@@H]([C@H](O)[C@H](O[C@@H]3O[C@@H]([C@@H](O)[C@H](O)[C@H]3NC(=O)C)CO)[C@H]2O)CO',
                                      'name': 'N-[(2S,3R,4R,5S,6R)-2-[(2S,3R,4S,5S,6R)-2-[(2R,3S,4R,5R,6R)-5-Acetamido-6-[(2R,3R,4R,5S)-5-acetamido-2,3,4,6-tetrahydroxyhexoxy]-4-hydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-4-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]acetamide',
                                      'reason': 'Found only 1 sugar rings, '
                                                'need at least 2'},
                                  {   'smiles': 'O=C(N[C@H]1[C@@H](O[C@H](CO)[C@H]([C@@H]1O)O[C@H]2O[C@@H]([C@@H](O)[C@@H]([C@H]2[N+](C)(C)C)O)CO)O[C@H]3[C@H](O)[C@@H](NC(=O)C)[C@H](O[C@H]4[C@H](O)C(NC(=O)C)C(O)O[C@@H]4CO)O[C@@H]3CO)C',
                                      'name': 'TMG-chitotriomycin',
                                      'reason': 'Found only 0 sugar rings, '
                                                'need at least 2'},
                                  {   'smiles': 'O=C1C2=C(O)C3=C(O)C(=C(O[C@@H]4O[C@@H]([C@H](OC(=O)C)[C@@H](C4)O[C@@H]5O[C@@H]([C@H](OC)[C@@H](C5)O)C)C)C=C3C=C2CC(C1O[C@@H]6O[C@@H]([C@@H](O)[C@@H](C6)O[C@@H]7O[C@@H]([C@@H](O)[C@@H](C7)O[C@@H]8O[C@H]([C@H](OC(=O)C)[C@](C8)(O)C)C)C)C)C(OC)C(=O)C(O)C(O)C)C',
                                      'name': '[(2R,3S,4R,6S)-6-[[6-[(2S,4R,5R,6R)-4-[(2S,4R,5R,6R)-4-[(2S,4S,5S,6S)-5-acetyloxy-4-hydroxy-4,6-dimethyloxan-2-yl]oxy-5-hydroxy-6-methyloxan-2-yl]oxy-5-hydroxy-6-methyloxan-2-yl]oxy-7-(3,4-dihydroxy-1-methoxy-2-oxopentyl)-4,10-dihydroxy-3-methyl-5-oxo-7,8-dihydro-6H-anthracen-2-yl]oxy]-4-[(2S,4R,5R,6R)-4-hydroxy-5-methoxy-6-methyloxan-2-yl]oxy-2-methyloxan-3-yl] '
                                              'acetate',
                                      'reason': 'Found only 0 sugar rings, '
                                                'need at least 2'},
                                  {   'smiles': 'O1[C@@H](OC(C[C@H](O)CCCC[C@@H](O)CC(O)=O)C)C(O)C[C@@H](O)C1C',
                                      'name': '(3R,8R)-10-[(2R,5R)-3,5-Dihydroxy-6-methyloxan-2-yl]oxy-3,8-dihydroxyundecanoic '
                                              'acid',
                                      'reason': 'Found only 0 sugar rings, '
                                                'need at least 2'},
                                  {   'smiles': 'O([C@H]1[C@@H](O)[C@H](O[C@@H](O)[C@@H]1O)CO[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2NC(=O)C)CO)[C@@H]3O[C@@H]([C@@H](O)[C@H](O)[C@H]3NC(=O)C)CO',
                                      'name': 'N-[(2R,3R,4R,5S,6R)-2-[[(2R,3S,4S,5R,6R)-4-[(2S,3R,4R,5S,6R)-3-Acetamido-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-3,5,6-trihydroxyoxan-2-yl]methoxy]-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]acetamide',
                                      'reason': 'Found only 1 sugar rings, '
                                                'need at least 2'},
                                  {   'smiles': 'O=C1C2=C(O)C3=C(O)C(=C(O[C@@H]4O[C@@H]([C@@H](O)[C@@H](C4)O)C)C=C3C=C2C[C@H]([C@@H]1O[C@@H]5O[C@@H]([C@@H](O)[C@@H](C5)O[C@@H]6O[C@@H]([C@H](O)[C@@H](C6)O[C@@H]7O[C@@H]([C@@H](O)[C@](C7)(O)C)C)C)C)[C@H](OC)C(=O)[C@@H](O)[C@H](O)C)C',
                                      'name': '3A-deolivosyl-mithramycin',
                                      'reason': 'Found only 0 sugar rings, '
                                                'need at least 2'}],
    'attempt': 1,
    'reasoning': 'Looking at the previous results, I notice several issues '
                 'that need to be addressed:\n'
                 '\n'
                 '1. The hydroxyl group detection was too strict - many valid '
                 'oligosaccharides were missed because the hydroxyl pattern '
                 "matching wasn't correctly identifying all OH groups. The "
                 'SMARTS pattern "[O][H]" may not work well since hydrogens '
                 'are often implicit in SMILES.\n'
                 '\n'
                 '2. The ring detection was too simplistic - the SMARTS '
                 "patterns for oxane and oxolane rings weren't specific enough "
                 'for sugar rings, which typically have specific substitution '
                 'patterns.\n'
                 '\n'
                 '3. Some valid oligosaccharides with modifications (like '
                 'acetyl groups) were being missed.\n'
                 '\n'
                 '4. The molecular weight thresholds may have been too '
                 'restrictive.\n'
                 '\n'
                 "Here's my improved approach:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 1428,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 100,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 0.9345549738219895,
    'f1': 0.9661705006765899,
    'accuracy': 0.9345549738219895,
    'negative_predictive_value': 0.0}