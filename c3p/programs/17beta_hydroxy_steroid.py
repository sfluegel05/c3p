"""
Classifies: CHEBI:35343 17beta-hydroxy steroid
"""
"""
Classifies: 17beta-hydroxy steroid
A 17-hydroxy steroid in which the hydroxy group at position 17 has a beta-configuration.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_17beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17beta-hydroxy steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 17beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for basic steroid core (four fused rings)
    steroid_core = Chem.MolFromSmarts("[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~1~[#6]~2~[#6]~3~4")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Check for 17-OH group in beta configuration
    # [C] is carbon 17, [OH1] is hydroxy group, '@' indicates stereochemistry
    # The [C] must be connected to 4 atoms (saturated)
    # Note: The exact SMARTS pattern depends on the numbering convention used
    oh_17_beta = Chem.MolFromSmarts('[C;X4](@[*])(@[*])(@[*])[OH1]')
    
    if not mol.HasSubstructMatch(oh_17_beta):
        return False, "No hydroxyl group with correct connectivity found"

    # Get matches for OH group
    oh_matches = mol.GetSubstructMatches(oh_17_beta)
    
    # Check if any of the matches are at position 17
    found_17_beta_oh = False
    for match in oh_matches:
        c_atom = mol.GetAtomWithIdx(match[0])  # Get the carbon atom
        # Check if this carbon is part of the D ring (ring 4) of the steroid
        # by checking its environment
        ring_info = mol.GetRingInfo()
        if ring_info.NumAtomRings(match[0]) > 0:  # Carbon must be part of a ring
            # Check chirality of the carbon
            if c_atom.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
                found_17_beta_oh = True
                break

    if not found_17_beta_oh:
        return False, "No 17-beta hydroxyl group found"

    # Additional validation: molecule should have reasonable size for a steroid
    num_atoms = mol.GetNumAtoms()
    if num_atoms < 20 or num_atoms > 100:
        return False, "Molecule size not consistent with steroid structure"

    # Count rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings for steroid structure"

    return True, "Contains steroid core with 17-beta hydroxyl group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35343',
                          'name': '17beta-hydroxy steroid',
                          'definition': 'A 17-hydroxy steroid in which the '
                                        'hydroxy group at position 17 has a '
                                        'beta-configuration.',
                          'parents': ['CHEBI:36838'],
                          'xrefs': [   'KEGG:C03051',
                                       'MetaCyc:17-beta-Hydroxysteroids'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'Cc1cc(C)c2Cc3c(C)cc(C)c(Cc4c(C)cc(C)c(Cc5c(C)cc(C)c(Cc1c2C)c5C)c4C)c3C',
                                     'name': '4,6,10,12,16,18,22,24,25,26,27,28-dodecamethylcalix[4]arene',
                                     'reason': 'No steroid core structure '
                                               'found'},
                                 {   'smiles': 'O=C1C(=C2C=C3[C@]([C@@H](C(C)C)[C@@H]([C@H]3O)OC(=O)C)(C)CC[C@]2(C)CC1)COC(=O)C',
                                     'name': 'Dahliane E',
                                     'reason': 'No steroid core structure '
                                               'found'},
                                 {   'smiles': 'O[C@@H]([C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)N[C@@H](CC(C)C)C(O)=O)C',
                                     'name': 'Glu-Thr-Leu',
                                     'reason': 'No steroid core structure '
                                               'found'},
                                 {   'smiles': 'CCOc1ccc(NC(=O)C(C)O)cc1',
                                     'name': 'p-Lactophenetide',
                                     'reason': 'No steroid core structure '
                                               'found'},
                                 {   'smiles': 'O=C1NCC=2C1=C3C(N([C@@H]4O[C@H]([C@@H](O)[C@H]([C@H]4OC)O)C)C5=C3C=CC=C5)=C6NC7=C(C26)C=CC=C7',
                                     'name': "3'-epi-5'-methoxy-K252d",
                                     'reason': 'No steroid core structure '
                                               'found'},
                                 {   'smiles': 'O=C(OC1C(O)C(OC(C1O)CO)OCC2OC(OCC(OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)C(O)C(C2O)O)CCCCCCCCCCCCCCC',
                                     'name': '[(2S)-2-hexadecanoyloxy-3-[(2S,3R,4S,5S,6R)-6-[[(2S,3R,4S,5S,6R)-4-hexadecanoyloxy-3,5-dihydroxy-6-(hydroxymethyl)tetrahydropyran-2-yl]oxymethyl]-3,4,5-trihydroxy-tetrahydropyran-2-yl]oxy-propyl] '
                                             '(9E,12E,15E)-octadeca-9,12,15-trienoate',
                                     'reason': 'No steroid core structure '
                                               'found'},
                                 {   'smiles': 'O=C1C2=C(O)C=C(OC)C=C2C(=O)C3=C1[C@@H]([C@@H](O)[C@]([C@@H]3O)(O)C)C',
                                     'name': 'Altersolanol G',
                                     'reason': 'No steroid core structure '
                                               'found'},
                                 {   'smiles': '[H][C@]1(O[C@](O)(C[C@H](O)[C@H]1NC(=O)CO)C([O-])=O)[C@H](O)[C@H](O)CO',
                                     'name': 'N-glycoloyl-alpha-neuraminate',
                                     'reason': 'No steroid core structure '
                                               'found'},
                                 {   'smiles': 'OC(C(O)C/C=C\\C/C=C\\C/C=C\\CC)C/C=C\\C/C=C\\CCC(O)=O',
                                     'name': '10,11-DiHDPE',
                                     'reason': 'No steroid core structure '
                                               'found'},
                                 {   'smiles': '[Na+].[H][C@]12SCC(C)=C(N1C(=O)[C@H]2NC(=O)[C@H](N)c1ccccc1)C([O-])=O',
                                     'name': 'cephalexin sodium',
                                     'reason': 'No steroid core structure '
                                               'found'}],
    'sample_false_negatives': [   {   'smiles': '[H][C@]12CCC(=O)C=C1CC[C@@]1([H])[C@]3([H])CC[C@@](O)(C#C)[C@@]3(CC)CC(=C)[C@]21[H]',
                                      'name': 'etonogestrel',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': 'C1=C(OS(O)(=O)=O)C=CC2=C1CC[C@]3([C@@]4(CC[C@]([C@]4(CC[C@@]32[H])C)(C#C)O)[H])[H]',
                                      'name': '17alpha-ethynylestradiol '
                                              '3-sulfate',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': '[H][C@@]12C[C@H](O)[C@H](O)[C@@]1(C)CC[C@@]1([H])[C@@]2([H])CC=C2C[C@H](O)CC[C@]12C',
                                      'name': 'androst-5-ene-3alpha,16beta,17beta-triol',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': 'C[C@@]12[C@]([C@]3([C@](CC1)([C@@]4(C(CC3)=C[C@H](CC4)O)[H])[H])[H])(CC[C@@H]2O)[H]',
                                      'name': 'bolandiol',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': 'C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)[C@H](F)C[C@]34C)[C@@H]1CC[C@@H]2O',
                                      'name': '2alpha-fluoro-17beta-hydroxyandrost-4-en-3-one',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': '[H][C@@]12CC[C@@]3([H])[C@]4([H])CC[C@H](O)[C@@]4(C)CC[C@]3([H])[C@@]1(C)CCC(=O)C2',
                                      'name': '17beta-hydroxy-5alpha-androstan-3-one',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': '[H][C@]12CC[C@]3(C)[C@@H](O)CC[C@@]3([H])[C@]1([H])CCc1cc(OC(=O)N(CCCl)CCCl)ccc21',
                                      'name': 'estramustine',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': '[H][C@@]12CCC3=CC(=O)CCC3=C1C=C[C@@]1(C)[C@@]2([H])CC[C@]1(C)O',
                                      'name': '17beta-hydroxy-17-methylestra-4,9,11-trien-3-one',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': 'C1C2=CC(C[C@@H]([C@@]2([C@]3(CC[C@]4([C@]([C@@]3(C1)[H])(CC[C@@H]4O)[H])C)[H])C)O)=O',
                                      'name': '1alpha-hydroxytestosterone',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': 'C[C@]12CC[C@H]3[C@@H](CC[C@H]4C[C@@H](CC[C@]34C)O[C@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)[C@@H]1CC[C@@H]2O',
                                      'name': '5alpha-androstane-3alpha-glcosyl-17beta-diol',
                                      'reason': 'No steroid core structure '
                                                'found'}],
    'attempt': 0,
    'reasoning': "I'll create a program to classify 17beta-hydroxy steroids. "
                 'The key characteristics are:\n'
                 '1. Must have a steroid core structure\n'
                 '2. Must have a hydroxy group at position 17\n'
                 '3. The hydroxy group must be in beta configuration\n'
                 '\n'
                 "Here's the implementation:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 142219,
    'num_false_negatives': 81,
    'num_negatives': None,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': 0.9994307800421645,
    'negative_predictive_value': 0.9994307800421645}