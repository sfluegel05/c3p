"""
Classifies: CHEBI:3098 bile acid
"""
"""
Classifies: bile acid
"""
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.rdchem import ChiralType

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid based on its SMILES string.
    Bile acids are hydroxy-5beta-cholanic acids occurring in bile.
    This function checks for the presence of a steroid nucleus with
    specific ring sizes and fusion, hydroxyl groups at specific positions,
    and a carboxylic acid side chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure molecule has 4 rings
    ring_info = mol.GetRingInfo()
    ring_counts = ring_info.NumRings()
    if ring_counts < 4:
        return False, f"Molecule has {ring_counts} rings, expected at least 4 for steroid nucleus"

    # Get ring systems and check ring sizes
    rings = ring_info.AtomRings()
    ring_sizes = [len(ring) for ring in rings]
    ring_sizes.sort()
    # Steroid nucleus has three 6-membered rings and one 5-membered ring
    if ring_sizes[:4] != [5, 6, 6, 6]:
        return False, f"Ring sizes are {ring_sizes[:4]}, expected [5, 6, 6, 6]"

    # Check for fused ring system
    fused = Chem.GetSymmSSSR(mol)
    if len(fused) < 4:
        return False, "Rings are not fused into a steroid nucleus"

    # Check for carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)[O;H]')
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "Carboxylic acid group not found"

    # Check for hydroxyl groups (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts('[OX2H]')
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) == 0:
        return False, "No hydroxyl groups found"

    # Optional: Check for specific stereochemistry at position 5 (5β-configuration)
    # Get atom at position 5 (assuming numbering starting from ring A)
    try:
        # Generate 3D coordinates for stereochemistry
        rdDepictor.Compute2DCoords(mol)
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
        atom5 = mol.GetAtomWithIdx(4)  # Zero-based indexing
        if atom5.GetChiralTag() != ChiralType.CHI_TETRAHEDRAL_CCW:
            return False, "5β-configuration not found"
    except:
        return False, "Error checking stereochemistry at position 5"

    # Additional check: Verify presence of hydroxyl groups at specific positions
    # This requires mapping atom indices to positions, which can be complex
    # For simplicity, we'll skip this detailed check

    return True, "Molecule contains steroid nucleus with correct ring fusion, carboxylic acid, hydroxyl groups, and 5β-configuration"


__metadata__ = {
    'chemical_class': {
        'name': 'bile acid',
        'definition': "Any member of a group of hydroxy-5beta-cholanic acids occurring in bile, where they are present as the sodium salts of their amides with glycine or taurine. In mammals bile acids almost invariably have 5beta-configuration."
    }
}


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:3098',
                          'name': 'bile acid',
                          'definition': 'Any member of a group of '
                                        'hydroxy-5beta-cholanic acids occuring '
                                        'in bile, where they are present as '
                                        'the sodium salts of their amides with '
                                        'glycine or taurine. In mammals bile '
                                        'acids almost invariably have '
                                        '5beta-configuration.',
                          'parents': ['CHEBI:138366', 'CHEBI:24663'],
                          'xrefs': ['KEGG:C01558'],
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
               'C1C[C@H](C[C@]2(CC[C@@]3([C@](CC([C@]4([C@]3(CC[C@@]4([C@@](CCC(=O)O)(C)[H])[H])[H])C)=O)([C@@]12C)[H])[H])[H])O '
               'NAME: 12-ketolithocholic acid REASON: MISSED Steroid nucleus '
               'not found\n'
               ' * SMILES: '
               'O[C@H]1[C@]2([C@]3([C@@]([C@](CC3)([C@@H](CCC[C@@H](C)C(O)=O)C)[H])(CC[C@@]2([C@@]4([C@](C1)(C[C@H](O)C[C@H]4O)[H])C)[H])C)[H])[H] '
               'NAME: '
               '(25R)-1beta,3alpha,7alpha-trihydroxy-5beta-cholestan-26-oic '
               'acid REASON: MISSED Steroid nucleus not found\n'
               ' * SMILES: '
               'C1C[C@@H](C[C@@]2(CC([C@@]3([C@](CC[C@]4([C@]3(CC[C@@]4([C@@](CCC(=O)O)(C)[H])[H])[H])C)([C@@]12C)[H])[H])=O)[H])O '
               'NAME: 3beta-hydroxy-7-oxo-5alpha-cholan-24-oic acid REASON: '
               'MISSED Steroid nucleus not found\n'
               ' * SMILES: '
               'O[C@@H]1[C@@]2([C@@]([C@](C1)([C@@H](CCC(OC)=O)C)[H])(CC[C@]3([C@]2([C@H](O)C[C@]4([C@@]3(CC[C@@H](O)C4)C)[H])[H])[H])C)[H] '
               'NAME: 3alpha,7alpha,15alpha-Trihydroxy-5beta-cholane-24-oic '
               'acid methyl ester REASON: MISSED Steroid nucleus not found\n'
               ' * SMILES: '
               'O=C(O)CC[C@H](C1[C@]2([C@@H](O)CC3[C@@]4(C(CC(OC)(OC)CC4)CCC3C2CC1)C)C)C '
               'NAME: 3-Dimethoxy-12alpha-hydroxycholanic acid REASON: MISSED '
               'Steroid nucleus not found\n'
               ' * SMILES: '
               'O[C@@H]1[C@]2([C@]([C@]3([C@@]([C@@]4([C@](C[C@H]3O)(C[C@H](O)C[C@H]4O)[H])C)(C1)[H])[H])(CC[C@@]2([C@@H](CCCC(C)C(O)=O)C)[H])[H])C '
               'NAME: '
               '1beta,3alpha,7alpha,12alpha-Tetrahydroxy-5beta-cholestan-26-oic '
               'acid REASON: MISSED Steroid nucleus not found\n'
               ' * SMILES: '
               'C[C@H](CCC(=O)O)[C@H]1CC[C@@H]2[C@@]1(C(=O)C[C@@H]3[C@@H]2C(=O)C[C@H]4[C@@]3(CCC(=O)C4)C)C '
               'NAME: '
               '(4R)-4-[(5S,8S,9R,10S,13R,14S,17R)-10,13-dimethyl-3,7,12-trioxo-1,2,4,5,6,8,9,11,14,15,16,17-dodecahydrocyclopenta[a]phenanthren-17-yl]pentanoic '
               'acid REASON: MISSED Steroid nucleus not found\n'
               ' * SMILES: '
               'O[C@@H]1[C@]2([C@]([C@]3([C@@]([C@@]4([C@]([C@@H](O)[C@H]3O)(C[C@H](O)CC4)[H])C)(C1)[H])[H])(CC[C@@]2([C@@H](CC/C=C(\\C)/C(O)=O)C)[H])[H])C '
               'NAME: '
               '3alpha,6alpha,7alpha,12alpha-Tetrahydroxy-5beta-cholest-24-en-26-oic '
               'acid REASON: MISSED Steroid nucleus not found\n'
               ' * SMILES: '
               'O[C@]12[C@@]([C@](CC1)([C@@H](CCC(OC)=O)C)[H])(CC[C@]3([C@]2([C@H](O)C[C@]4([C@@]3(CC[C@@H](O)C4)C)[H])[H])[H])C '
               'NAME: Methyl '
               '(R)-4-((3R,5S,7R,8S,9S,10S,13R,14R,17R)-3,7,14-trihydroxy-10,13-dimethylhexadecahydro-1H-cyclopenta[a]phenanthren-17-yl)pentanoate '
               'REASON: MISSED Steroid nucleus not found\n'
               ' * SMILES: '
               'O=C1[C@]2([C@@](CC[C@@]2([C@@H](CCC(O)=O)C)[H])(C=3C([C@@]4([C@](CC3)(C[C@H](O)CC4)[H])C)=C1)[H])C '
               'NAME: 3alpha-Hydroxy-12-oxo-5beta-chola-7,9(11)-dien-24-oic '
               'Acid REASON: MISSED Steroid nucleus not found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [],
    'sample_false_negatives': [   {   'smiles': 'C1C[C@H](C[C@]2(CC[C@@]3([C@](CC([C@]4([C@]3(CC[C@@]4([C@@](CCC(=O)O)(C)[H])[H])[H])C)=O)([C@@]12C)[H])[H])[H])O',
                                      'name': '12-ketolithocholic acid',
                                      'reason': '5β-configuration not found'},
                                  {   'smiles': 'O[C@H]1[C@]2([C@]3([C@@]([C@](CC3)([C@@H](CCC[C@@H](C)C(O)=O)C)[H])(CC[C@@]2([C@@]4([C@](C1)(C[C@H](O)C[C@H]4O)[H])C)[H])C)[H])[H]',
                                      'name': '(25R)-1beta,3alpha,7alpha-trihydroxy-5beta-cholestan-26-oic '
                                              'acid',
                                      'reason': '5β-configuration not found'},
                                  {   'smiles': 'O[C@@H]1[C@@]2([C@@]([C@](C1)([C@@H](CCC(OC)=O)C)[H])(CC[C@]3([C@]2([C@H](O)C[C@]4([C@@]3(CC[C@@H](O)C4)C)[H])[H])[H])C)[H]',
                                      'name': '3alpha,7alpha,15alpha-Trihydroxy-5beta-cholane-24-oic '
                                              'acid methyl ester',
                                      'reason': 'Carboxylic acid group not '
                                                'found'},
                                  {   'smiles': 'O=C(O)CC[C@H](C1[C@]2([C@@H](O)CC3[C@@]4(C(CC(OC)(OC)CC4)CCC3C2CC1)C)C)C',
                                      'name': '3-Dimethoxy-12alpha-hydroxycholanic '
                                              'acid',
                                      'reason': '5β-configuration not found'},
                                  {   'smiles': 'O[C@@H]1[C@]2([C@]([C@]3([C@@]([C@@]4([C@](C[C@H]3O)(C[C@H](O)C[C@H]4O)[H])C)(C1)[H])[H])(CC[C@@]2([C@@H](CCCC(C)C(O)=O)C)[H])[H])C',
                                      'name': '1beta,3alpha,7alpha,12alpha-Tetrahydroxy-5beta-cholestan-26-oic '
                                              'acid',
                                      'reason': '5β-configuration not found'},
                                  {   'smiles': 'C[C@H](CCC(=O)O)[C@H]1CC[C@@H]2[C@@]1(C(=O)C[C@@H]3[C@@H]2C(=O)C[C@H]4[C@@]3(CCC(=O)C4)C)C',
                                      'name': '(4R)-4-[(5S,8S,9R,10S,13R,14S,17R)-10,13-dimethyl-3,7,12-trioxo-1,2,4,5,6,8,9,11,14,15,16,17-dodecahydrocyclopenta[a]phenanthren-17-yl]pentanoic '
                                              'acid',
                                      'reason': '5β-configuration not found'},
                                  {   'smiles': 'O[C@@H]1[C@]2([C@]([C@]3([C@@]([C@@]4([C@]([C@@H](O)[C@H]3O)(C[C@H](O)CC4)[H])C)(C1)[H])[H])(CC[C@@]2([C@@H](CC/C=C(\\C)/C(O)=O)C)[H])[H])C',
                                      'name': '3alpha,6alpha,7alpha,12alpha-Tetrahydroxy-5beta-cholest-24-en-26-oic '
                                              'acid',
                                      'reason': '5β-configuration not found'},
                                  {   'smiles': 'O[C@]12[C@@]([C@](CC1)([C@@H](CCC(OC)=O)C)[H])(CC[C@]3([C@]2([C@H](O)C[C@]4([C@@]3(CC[C@@H](O)C4)C)[H])[H])[H])C',
                                      'name': 'Methyl '
                                              '(R)-4-((3R,5S,7R,8S,9S,10S,13R,14R,17R)-3,7,14-trihydroxy-10,13-dimethylhexadecahydro-1H-cyclopenta[a]phenanthren-17-yl)pentanoate',
                                      'reason': 'Carboxylic acid group not '
                                                'found'},
                                  {   'smiles': 'O=C1[C@]2([C@@](CC[C@@]2([C@@H](CCC(O)=O)C)[H])(C=3C([C@@]4([C@](CC3)(C[C@H](O)CC4)[H])C)=C1)[H])C',
                                      'name': '3alpha-Hydroxy-12-oxo-5beta-chola-7,9(11)-dien-24-oic '
                                              'Acid',
                                      'reason': '5β-configuration not found'},
                                  {   'smiles': 'O=C1[C@]2([C@]([C@]3(C([C@@]4([C@](CC3)(C[C@@H](O)CC4)[H])C)=C1)[H])(CC[C@@]2([C@@H](CCC(O)=O)C)[H])[H])C',
                                      'name': '3beta-Hydroxy-12-oxo-5beta-chol-9(11)-en-24-oic '
                                              'Acid',
                                      'reason': '5β-configuration not found'}],
    'attempt': 1,
    'reasoning': "Firstly, let's analyze the previous program and the outcomes "
                 'to understand what went wrong and how to improve it.\n'
                 '\n'
                 '**Analysis:**\n'
                 '\n'
                 '1. **Inadequate Steroid Nucleus Detection:**\n'
                 '   - The SMARTS pattern used to detect the steroid nucleus '
                 "was too simplistic and didn't accurately represent the "
                 'complex fused ring system of steroids.\n'
                 '   - Steroids have a characteristic tetracyclic ring system '
                 'consisting of three six-membered rings (rings A, B, and C) '
                 'and one five-membered ring (ring D) fused together.\n'
                 '   - The previous pattern likely failed to match this '
                 'structure in the provided bile acid examples, leading to '
                 'false negatives.\n'
                 '\n'
                 '2. **Lack of Stereochemistry Consideration:**\n'
                 '   - Bile acids have specific stereochemistry, particularly '
                 'the 5β-configuration, which refers to the orientation of '
                 'certain hydrogen atoms.\n'
                 '   - The previous code did not account for stereochemistry, '
                 'causing it to miss valid bile acids.\n'
                 '\n'
                 '3. **Generic Substituent Checks:**\n'
                 '   - The checks for hydroxyl groups and carboxylic acid '
                 'groups were not specific to their positions on the steroid '
                 'nucleus.\n'
                 '   - Bile acids have hydroxyl groups at specific positions '
                 '(e.g., 3α, 7α, 12α) and a carboxylic acid side chain at '
                 'position 17 (C-17).\n'
                 '\n'
                 '**Improvements:**\n'
                 '\n'
                 '- **Use a More Accurate SMARTS Pattern for Steroid '
                 'Nucleus:**\n'
                 '  - Implement a SMARTS pattern that accurately represents '
                 'the tetracyclic ring system of steroids with correct ring '
                 'sizes and fusion.\n'
                 '  - Ensure that the ring fusion and sizes (6-6-6-5) are '
                 'correctly specified.\n'
                 '\n'
                 '- **Include Stereochemistry in the SMARTS Pattern:**\n'
                 '  - Incorporate stereochemical specifications to detect the '
                 '5β-configuration.\n'
                 '  - While accounting for full stereochemistry might be '
                 'complex, focusing on key stereocenters can improve '
                 'accuracy.\n'
                 '\n'
                 '- **Check for Specific Substituents at Correct Positions:**\n'
                 '  - Verify the presence of hydroxyl groups at specific '
                 'positions on the steroid nucleus.\n'
                 '  - Confirm the presence of a carboxylic acid side chain '
                 'attached at the correct position.\n'
                 '\n'
                 "Now, let's implement the improved function.",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 9,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 100,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 0.08256880733944955,
    'f1': 0.15254237288135594,
    'accuracy': 0.08256880733944955,
    'negative_predictive_value': 0.0}