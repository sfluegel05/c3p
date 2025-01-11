"""
Classifies: CHEBI:50126 tetrasaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tetrasaccharide(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is a tetrasaccharide based on its SMILES string.
    A tetrasaccharide consists of exactly 4 monosaccharide units connected by glycosidic bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple[bool, str]: (is_tetrasaccharide, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic checks for carbohydrate characteristics
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 20:  # Typical tetrasaccharide has at least 20 carbons
        return False, "Too few carbons for a tetrasaccharide"
    if o_count < 15:  # Typical tetrasaccharide has at least 15 oxygens
        return False, "Too few oxygens for a tetrasaccharide"

    # More flexible pattern for sugar rings (both pyranose and furanose)
    pyranose_pattern = Chem.MolFromSmarts("[C]1[C][C][C]([O,C])[C]O1")
    furanose_pattern = Chem.MolFromSmarts("[C]1[C][C]([O,C])[C]O1")
    
    pyranose_matches = len(mol.GetSubstructMatches(pyranose_pattern))
    furanose_matches = len(mol.GetSubstructMatches(furanose_pattern))
    total_sugar_rings = pyranose_matches + furanose_matches
    
    if total_sugar_rings < 4:
        return False, f"Found only {total_sugar_rings} sugar rings, need at least 4"
    
    # Look for glycosidic bonds between sugar rings
    # Pattern matches C-O-C where one C is part of a ring
    glycosidic_pattern = Chem.MolFromSmarts("[C;R][O][C]")
    glycosidic_matches = len(mol.GetSubstructMatches(glycosidic_pattern))
    
    if glycosidic_matches < 3:
        return False, "Too few glycosidic bonds for a tetrasaccharide"

    # Count oxygen atoms bonded to carbon (better than counting explicit hydroxyls)
    c_o_pattern = Chem.MolFromSmarts("[C][O]")
    c_o_count = len(mol.GetSubstructMatches(c_o_pattern))
    
    if c_o_count < 15:  # Tetrasaccharides typically have many C-O bonds
        return False, "Too few C-O bonds for a tetrasaccharide"

    # Check for characteristics that would exclude it from being a tetrasaccharide
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    if aromatic_atoms > 2:  # Allow some aromatic atoms for modified sugars
        return False, "Too many aromatic atoms for a tetrasaccharide"

    # Look for typical sugar carbon patterns
    sugar_carbon_pattern = Chem.MolFromSmarts("[C]([O])([O])[C]([O])")
    if not mol.HasSubstructMatch(sugar_carbon_pattern):
        return False, "Missing typical sugar carbon pattern"

    # Additional check for anomeric carbons (carbons with two oxygens attached)
    anomeric_pattern = Chem.MolFromSmarts("[C;R]([O])([O])")
    anomeric_matches = len(mol.GetSubstructMatches(anomeric_pattern))
    
    if anomeric_matches < 3:  # Should have at least 3 anomeric carbons
        return False, "Too few anomeric carbons for a tetrasaccharide"

    return True, "Contains 4 monosaccharide units connected by glycosidic bonds"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:50126',
                          'name': 'tetrasaccharide',
                          'definition': 'An oligosaccharide comprising four '
                                        'monomeric monosaccharide units.',
                          'parents': ['CHEBI:50699'],
                          'xrefs': ['KEGG:C06768'],
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
               'C[C@@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@@H](O)[C@@H](CO)O[C@H]2O[C@@H]2[C@@H](CO)OC(O)[C@H](O)[C@H]2O[C@@H]2O[C@@H](C)[C@@H](O)[C@@H](O)[C@@H]2O)[C@@H](O)[C@H](O)[C@@H]1O '
               'NAME: '
               'alpha-L-Fucp-(1->2)-beta-D-Galp-(1->4)-[alpha-L-Fucp-(1->3)]-D-Glcp '
               'REASON: MISSED Too few hydroxyl groups for a tetrasaccharide\n'
               ' * SMILES: '
               'OC[C@H]1O[C@H](O[C@H]2[C@@H](O)[C@H](O)[C@@H](CO)O[C@@H]2OC[C@H]2O[C@H](OC[C@H]3OC(O)[C@@H](O)[C@@H](O)[C@@H]3O)[C@@H](O)[C@@H](O)[C@@H]2O)[C@@H](O)[C@@H](O)[C@@H]1O '
               'NAME: '
               'alpha-D-Manp-(1->2)-alpha-D-Manp-(1->6)-alpha-D-Manp-(1->6)-D-Manp '
               'REASON: MISSED Too few hydroxyl groups for a tetrasaccharide\n'
               ' * SMILES: '
               'O[C@H]1[C@@H]([C@@H]([C@@H](O[C@@H]1CO[C@H]2[C@@H]([C@H]([C@H]([C@H](O2)CO)O)O)O)O[C@H]3[C@@H]([C@@H]([C@@H](O[C@@H]3CO)O[C@H]4[C@@H]([C@@H](C(O[C@@H]4CO)O)O)O)O)O)O)O '
               'NAME: '
               'alpha-D-Galp-(1->6)-beta-D-Manp-(1->4)-beta-D-Manp-(1->4)-D-Manp '
               'REASON: MISSED Too few hydroxyl groups for a tetrasaccharide\n'
               ' * SMILES: '
               'C[C@@H]1O[C@@H](O[C@H](CO)[C@@H](O)[C@@H](O)CO)[C@H](O)[C@H](O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O[C@H]3O[C@H](CO)[C@H](O)[C@H](O)[C@H]3O)[C@H]2O)[C@H]1O '
               'NAME: '
               'alpha-D-Galp-(1->3)-alpha-D-Glcp-(1->3)-alpha-L-Rhap-(1->4)-D-ribitol '
               'REASON: MISSED Too few hydroxyl groups for a tetrasaccharide\n'
               ' * SMILES: '
               'C[C@H]1O[C@H](O[C@H]2[C@H](O)[C@@H](CO)O[C@H](O[C@H]3[C@H](C)O[C@@H](O[C@@H]4[C@@H](O)[C@@H](O)O[C@H](CO)[C@@H]4O)[C@H](O)[C@@H]3O)[C@H]2O)[C@@H](O)C[C@@H]1O '
               'NAME: '
               'alpha-Tyvp-(1->3)-alpha-D-Manp-(1->4)-alpha-L-Rhap-(1->3)-alpha-D-Galp '
               'REASON: MISSED Too few hydroxyl groups for a tetrasaccharide\n'
               ' * SMILES: '
               'C[C@@H]1O[C@@H](O[C@H]2[C@@H](O[C@H](CO)[C@H](O)[C@@H]2O[C@H]2O[C@H](CO)[C@H](O)[C@H](O)[C@H]2O)O[C@@H]2[C@@H](CO)O[C@@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@H](O)[C@@H]1O '
               'NAME: '
               'alpha-D-Galp-(1->3)-[alpha-L-Fucp-(1->2)]-beta-D-Galp-(1->4)-beta-D-Glcp '
               'REASON: MISSED Too few hydroxyl groups for a tetrasaccharide\n'
               ' * SMILES: '
               '[H][C@]1(O[C@@H](O[C@H](CO)[C@]2([H])O[C@@H](O[C@H](CO)[C@]3([H])O[C@@H](O[C@H](CO)[C@]4([H])O[C@@H](O)[C@H](O)[C@H]4O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@H]1O)[C@H](O)CO '
               'NAME: '
               'beta-D-Galf-(1->5)-beta-D-Galf-(1->5)-beta-D-Galf-(1->5)-beta-D-Galf '
               'REASON: MISSED Found only 0 sugar rings, need exactly 4\n'
               ' * SMILES: O[C@@H]([C@@H](O)C(O)=O)C NAME: 4-Deoxyerythronic '
               'acid REASON: MISSED Too few carbons for a tetrasaccharide\n'
               ' * SMILES: '
               'O([C@H]1[C@H](O)[C@@H](O[C@@H]2O[C@H]([C@@H](O)[C@@H](O)[C@@H]2O)C)C(O[C@@H]1CO)O)[C@@H]3O[C@@H]([C@H](O)[C@H](O)[C@H]3O[C@@H]4O[C@H]([C@@H](O)[C@@H](O)[C@@H]4O)C)CO '
               "NAME: 2,2'-difucosyllactose REASON: MISSED Too few hydroxyl "
               'groups for a tetrasaccharide\n'
               ' * SMILES: '
               'OC[C@H]1O[C@H](O[C@H]2[C@H](O)[C@@H](O)[C@@H](O[C@H]3[C@H](O)[C@@H](O)[C@@H](O[C@H]4[C@H](O)[C@@H](O)[C@H](O)O[C@@H]4CO)O[C@@H]3CO)O[C@@H]2CO)[C@H](O)[C@@H](O)[C@@H]1O '
               'NAME: '
               'alpha-D-Glcp-(1->4)-alpha-D-Glcp-(1->4)-alpha-D-Glcp-(1->4)-beta-D-Glcp '
               'REASON: MISSED Too few hydroxyl groups for a tetrasaccharide\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C1C=C(N=C2N1NN=N2)C3=CC=CC=C3',
                                     'name': '5-phenyl-1,7-dihydrotetrazolo[1,5-a]pyrimidine',
                                     'reason': 'Too few carbons for a '
                                               'tetrasaccharide'},
                                 {   'smiles': 'O(C(=O)CCCCCCCCCCC/C=C\\CCCCCCCC)[C@H](COC(=O)CCCCCCCCC/C=C\\CCCCCC)COC(=O)CCCCCCC/C=C\\CCCC',
                                     'name': 'TG(14:1(9Z)/22:1(13Z)/18:1(11Z))',
                                     'reason': 'Too few oxygens for a '
                                               'tetrasaccharide'},
                                 {   'smiles': '[C@@]12([C@@]([C@]([C@@H](CC1)C)(CCC(CC)C)C)(CCC[C@@H]2C)[H])C',
                                     'name': 'clerodane',
                                     'reason': 'Too few oxygens for a '
                                               'tetrasaccharide'},
                                 {   'smiles': 'S(O[C@@H]1[C@H](O)[C@@H](O)[C@H](O[C@@H]([C@@H](O)[C@H](O)CO[C@]2(O[C@H]([C@H](NC(=O)C)[C@@H](O)C2)[C@H](O)[C@H](O)CO)C(O)=O)[C@@H](NC(=O)C)CO)O[C@@H]1CO)(O)(=O)=O',
                                     'name': '(2R,4S,5R,6R)-5-Acetamido-2-[(2R,3S,4R,5S)-5-acetamido-4-[(2R,3R,4R,5R,6R)-3,4-dihydroxy-6-(hydroxymethyl)-5-sulfooxyoxan-2-yl]oxy-2,3,6-trihydroxyhexoxy]-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'Found only 2 sugar rings, need '
                                               'at least 4'},
                                 {   'smiles': 'O=C(NC(CC(C)C)C(=O)N[C@@H](CCCN=C(N)N)C=O)C(NC(=O)CC)CC(C)C',
                                     'name': 'Leupeptin Pr-LL',
                                     'reason': 'Too few oxygens for a '
                                               'tetrasaccharide'},
                                 {   'smiles': 'COC1=CC=C(C=C1)CCNC(=O)C(=C2C3=CC=CC=C3C(=N2)NC(=O)C4=CC=CS4)C#N',
                                     'name': 'N-[3-[1-cyano-2-[2-(4-methoxyphenyl)ethylamino]-2-oxoethylidene]-1-isoindolyl]-2-thiophenecarboxamide',
                                     'reason': 'Too few oxygens for a '
                                               'tetrasaccharide'},
                                 {   'smiles': 'O=C1O[C@@H](CC=2C1=C(OC)C(O)=CC2)CCC[C@@H](O)C',
                                     'name': 'Penicimarin C',
                                     'reason': 'Too few carbons for a '
                                               'tetrasaccharide'},
                                 {   'smiles': 'CO[C@@H]1[C@@H]2[C@@H](C[C@@H]3[C@@](O2)([C@H]([C@@H]([C@H](O3)C(=O)OC)O)O)O)OC1N4C=NC5=C4N=CN=C5N',
                                     'name': 'LSM-4497',
                                     'reason': 'Too few carbons for a '
                                               'tetrasaccharide'},
                                 {   'smiles': 'C[C@@H]1[C@H](O)CCC2=CC[C@H](C[C@]12C)C(C)=C',
                                     'name': '1-deoxycapsidiol',
                                     'reason': 'Too few carbons for a '
                                               'tetrasaccharide'},
                                 {   'smiles': 'C1CN=C(Cc2ccccc2)N1',
                                     'name': 'tolazoline',
                                     'reason': 'Too few carbons for a '
                                               'tetrasaccharide'}],
    'sample_false_negatives': [   {   'smiles': 'O[C@@H]([C@@H](O)C(O)=O)C',
                                      'name': '4-Deoxyerythronic acid',
                                      'reason': 'Too few carbons for a '
                                                'tetrasaccharide'},
                                  {   'smiles': 'OCC(O)C(O)CC(O)=O',
                                      'name': '2-Deoxypentonic acid',
                                      'reason': 'Too few carbons for a '
                                                'tetrasaccharide'},
                                  {   'smiles': 'O(C(O[C@@H](C1=CC=CC=C1)C#N)C(O)C(OC)=O)C(CO)C(OC)=O',
                                      'name': 'Sambucus nigra Degraded '
                                              'cyanogenic glycosides '
                                              "(2'-Epimer)",
                                      'reason': 'Too few carbons for a '
                                                'tetrasaccharide'},
                                  {   'smiles': 'OC[C@@H](O)[C@H]([C@H](O)[C@@H](O)C(=O)[O-])O.[Zn+2].OC[C@@H](O)[C@H]([C@H](O)[C@@H](O)C(=O)[O-])O',
                                      'name': 'Zinc gluconate',
                                      'reason': 'Too few carbons for a '
                                                'tetrasaccharide'},
                                  {   'smiles': 'O[C@H]([C@H](O)CO)[C@@H](O)C([O-])=O.[Li+]',
                                      'name': 'D-ribonate lithium salt',
                                      'reason': 'Too few carbons for a '
                                                'tetrasaccharide'},
                                  {   'smiles': 'O([C@H]1[C@@H]([C@H]([C@H](O[C@H]2[C@@H]([C@H]([C@H](O[C@H]3[C@@H]([C@H]([C@H](O)[C@H](O3)CO)O)NC(=O)C)[C@H](O2)CO)O)NC(=O)C)[C@H](O1)CO)O)NC(=O)C)[C@@H]4[C@H](OC(=O)[C@@H]([C@H]4O)NC(=O)C)CO',
                                      'name': "N,N',N'',N'''-tetraacetylchitotetraono-1,5-lactone",
                                      'reason': 'Missing typical sugar carbon '
                                                'pattern'},
                                  {   'smiles': '[Ca+2].O[C@H]([C@@H](O)C([O-])=O)CO.O[C@H]([C@@H](O)C([O-])=O)CO',
                                      'name': 'L-Threonic acid Calcium Salt',
                                      'reason': 'Too few carbons for a '
                                                'tetrasaccharide'},
                                  {   'smiles': 'CC(C)=CCc1cc(cc(CC=C(C)C)c1O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)C(=O)OC[C@H]1CCN2CCC[C@H]12',
                                      'name': 'Auriculine',
                                      'reason': 'Too few oxygens for a '
                                                'tetrasaccharide'},
                                  {   'smiles': '[H][C@]1(O[C@H](COS(O)(=O)=O)[C@@H](O)[C@H](O)[C@H]1O)Oc1ccc(cc1)C(O)=O',
                                      'name': '4-hydroxybenzoic acid '
                                              '4-(6-O-sulfo)glucopyranoside',
                                      'reason': 'Too few carbons for a '
                                                'tetrasaccharide'},
                                  {   'smiles': 'O(C(=O)CCCCCCCC/C=C\\CCC)C[C@H](O)C(O)=O',
                                      'name': '2S-hydroxy-3-(10Z-tetradecenoyloxy)-propanoic '
                                              'acid',
                                      'reason': 'Too few carbons for a '
                                                'tetrasaccharide'}],
    'attempt': 1,
    'reasoning': 'Let me analyze the issues with the previous code and suggest '
                 'improvements:\n'
                 '\n'
                 '1. Ring Detection Problems:\n'
                 '- The previous sugar ring pattern '
                 '`[C]1[C][C][C]([O,C])[C]([O,C])O1` was too restrictive\n'
                 '- It missed furanose rings and some pyranose variations\n'
                 '- We should use a more flexible pattern that can detect both '
                 '5 and 6-membered sugar rings\n'
                 '\n'
                 '2. Hydroxyl Group Detection:\n'
                 '- The hydroxyl pattern `[O][H]` was problematic because '
                 "RDKit's SMILES parser often doesn't include explicit "
                 'hydrogens\n'
                 '- Many valid tetrasaccharides were rejected due to "Too few '
                 'hydroxyl groups"\n'
                 '- We should count oxygen atoms connected to carbons instead\n'
                 '\n'
                 '3. Ring Count Issues:\n'
                 "- Simply counting rings isn't sufficient as some sugars can "
                 'have additional small rings\n'
                 '- We need to specifically identify sugar rings rather than '
                 'just any rings\n'
                 '\n'
                 '4. Glycosidic Bond Detection:\n'
                 '- The pattern `[C][O][C]` was too general and could match '
                 'non-glycosidic C-O-C bonds\n'
                 '- We need a more specific pattern for glycosidic linkages '
                 'between sugar rings\n'
                 '\n'
                 "Here's an improved version:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 69,
    'num_false_positives': 100,
    'num_true_negatives': 1748,
    'num_false_negatives': 19,
    'num_negatives': None,
    'precision': 0.40828402366863903,
    'recall': 0.7840909090909091,
    'f1': 0.5369649805447471,
    'accuracy': 0.9385330578512396,
    'negative_predictive_value': 0.989247311827957}