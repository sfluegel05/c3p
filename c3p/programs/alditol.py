"""
Classifies: CHEBI:17522 alditol
"""
"""
Classifies: CHEBI:15972 alditol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_alditol(smiles: str):
    """
    Determines if a molecule is an alditol based on its SMILES string.
    An alditol is an acyclic polyol with formula HOCH2[CH(OH)]nCH2OH.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alditol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # First check for basic requirements
    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 3:
        return False, "Too few carbons for an alditol (minimum 3)"
    if o_count < 3:
        return False, "Too few oxygens for an alditol (minimum 3)"

    # Check for absence of carbonyl groups
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    if mol.HasSubstructMatch(carbonyl_pattern):
        return False, "Contains carbonyl group(s)"

    # Look for terminal CH2OH groups
    terminal_oh_pattern = Chem.MolFromSmarts("[CH2X4][OX2H1]")
    terminal_oh_matches = mol.GetSubstructMatches(terminal_oh_pattern)
    
    # For pure alditols, we need exactly two terminal CH2OH groups
    # However, some derivatives may have modifications, so we'll be lenient here
    if len(terminal_oh_matches) < 1:
        return False, "Missing terminal CH2OH group"

    # Check for chain of carbons with OH groups
    chain_with_oh_pattern = Chem.MolFromSmarts("[CH2X4,CHX4][OX2H1]")
    oh_chain_matches = mol.GetSubstructMatches(chain_with_oh_pattern)
    
    if len(oh_chain_matches) < 3:
        return False, "Insufficient hydroxyl groups in carbon chain"

    # For pure alditols, check if the molecule is mostly linear
    # by comparing the number of branches to the number of atoms
    n_atoms = mol.GetNumAtoms()
    n_branches = 0
    for atom in mol.GetAtoms():
        if atom.GetDegree() > 2:
            n_branches += 1
    
    # Calculate ratio of branching points to total atoms
    branch_ratio = n_branches / n_atoms if n_atoms > 0 else 1
    
    # If it's a pure alditol (not a derivative), it should be mostly linear
    if branch_ratio > 0.5:
        return False, "Too many branches for an alditol"

    # Check for presence of atoms other than C, H, O
    # (allowing for some derivatives that might contain other atoms)
    non_cho_atoms = sum(1 for atom in mol.GetAtoms() 
                       if atom.GetAtomicNum() not in [1, 6, 8])
    
    if non_cho_atoms > 0:
        # This might be an alditol derivative
        return True, "Appears to be an alditol derivative"
    
    # If we've made it here, check the ratio of OH groups to carbons
    # For pure alditols, this should be close to 1
    oh_ratio = o_count / c_count
    
    if 0.8 <= oh_ratio <= 1.2:
        return True, "Pure alditol structure detected"
    else:
        return True, "Possible alditol derivative"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:17522',
                          'name': 'alditol',
                          'definition': 'A carbohydrate that is an acyclic '
                                        'polyol having the general formula '
                                        'HOCH2[CH(OH)]nCH2OH (formally '
                                        'derivable from an aldose by reduction '
                                        'of the carbonyl group).',
                          'parents': ['CHEBI:16646', 'CHEBI:26191'],
                          'xrefs': ['KEGG:C00717', 'Wikipedia:Glycerin'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'C[C@@H](C1=CC=CC=C1)NC(=O)C[C@H]2CC[C@@H]([C@@H](O2)CO)NC(=O)CN3CCOCC3',
                                     'name': '2-[(2R,5S,6R)-6-(hydroxymethyl)-5-[[2-(4-morpholinyl)-1-oxoethyl]amino]-2-oxanyl]-N-[(1S)-1-phenylethyl]acetamide',
                                     'reason': 'Contains carbonyl group(s)'},
                                 {   'smiles': 'C[N+](C)(C)[C@@H](Cc1c[nH]c(n1)S(=O)C[C@H](NC(=O)CC[C@H]([NH3+])C([O-])=O)C([O-])=O)C([O-])=O',
                                     'name': 'N(alpha)-(L-gamma-glutamyl)-hercynyl-L-cysteine '
                                             'sulfoxide(1-)',
                                     'reason': 'Contains carbonyl group(s)'},
                                 {   'smiles': 'O(C1=CC=2[C@]3([C@](N(CC3)C)(N(C2C=C1)C)[H])C)C(=O)N4CCC=5C(C4)=CC=CC5',
                                     'name': 'quilostigmine',
                                     'reason': 'Too few oxygens for an alditol '
                                               '(minimum 3)'},
                                 {   'smiles': 'O[C@@H]1[C@]23[C@@]4(N(C[C@@]([C@]2(C[C@@]4([C@]56[C@]3(CC(=O)[C@](C5)(C([C@H]6O)=C)[H])[H])[H])[H])(CC1)C)CC)[H]',
                                     'name': 'Bullatine G',
                                     'reason': 'Contains carbonyl group(s)'},
                                 {   'smiles': 'NC1=NC=NC2=C1N=CN2[C@@H]3O[C@H](COP(=O)(O)O)[C@@H](OC(=O)[C@@H](N)CCC(O)=O)[C@H]3O',
                                     'name': "3'-L-glutamyl-AMP",
                                     'reason': 'Contains carbonyl group(s)'},
                                 {   'smiles': 'O1C2(C(C3(C(C4(C(CC3OC(=O)C)C(OC(=O)CC4)(C)C)C)CC2)C)CC15C6N(C=7C5=CC=CC7)C(=O)C(N6)C)C',
                                     'name': 'Teraspiridole C_130091',
                                     'reason': 'Contains carbonyl group(s)'},
                                 {   'smiles': 'O(C1[C@@H](OC(=O)C)C(O[C@@H](OC2=C(OC3=C(C2=O)C(O)=CC(O[C@@H]4OC([C@@H](O)[C@H](O)C4O)CO)=C3CC=C(C)C)C5=CC=C(OC)C=C5)[C@H]1O)C)[C@@H]6OC[C@@H](O)[C@H](OC(=O)C)C6O',
                                     'name': 'Sempervirenoside A',
                                     'reason': 'Contains carbonyl group(s)'},
                                 {   'smiles': 'COC[C@]1(C(=O)C2CCN1CC2)CO',
                                     'name': '(2S)-2-(hydroxymethyl)-2-(methoxymethyl)-1-azabicyclo[2.2.2]octan-3-one',
                                     'reason': 'Contains carbonyl group(s)'},
                                 {   'smiles': 'Oc1c(C2CC(Cc3ccccc23)c2ccc(OCc3ccc(cc3)C(F)(F)F)cc2)c(=O)oc2ccccc12',
                                     'name': 'Flocoumafen',
                                     'reason': 'Missing terminal CH2OH group'},
                                 {   'smiles': 'O[C@H]1CC=2C(N(C=3C1=CC=CC3)C(=O)N)=CC=CC2',
                                     'name': '(S)-MHD',
                                     'reason': 'Too few oxygens for an alditol '
                                               '(minimum 3)'}],
    'sample_false_negatives': [   {   'smiles': 'CC(=O)N[C@@H]1[C@@H](O)[C@H](O[C@@H]2O[C@H](CO)[C@H](O)[C@H](O)[C@H]2O)[C@@H](CO)O[C@H]1OC[C@@H](O)[C@H](O)[C@H](O)[C@@H](O)CO',
                                      'name': 'beta-D-Gal-(1->4)-beta-D-GlcNAc-(1->6)-D-Gal-OH',
                                      'reason': 'Contains carbonyl group(s)'},
                                  {   'smiles': 'O([C@H]1[C@H](O[C@@H]2O[C@H]([C@@H](O)[C@@H](O)[C@@H]2O)C)[C@H](O[C@@H](O[C@H]3[C@@H](O)[C@H](O[C@@H](OCC(O)CO)[C@@H]3O)CO)[C@@H]1NC(=O)C)CO)[C@@H]4O[C@@H]([C@H](O)[C@H](O[C@H]5O[C@@H]([C@H](O)[C@H](O)[C@H]5NC(=O)C)CO)[C@H]4O[C@@H]6O[C@H]([C@@H](O)[C@@H](O)[C@@H]6O)C)CO',
                                      'name': 'N-[(2R,3R,4R,5R,6R)-2-[(2R,3R,4S,5S,6R)-2-[(2S,3R,4R,5S,6R)-3-Acetamido-2-[(2R,3R,4S,5S,6R)-2-(2,3-dihydroxypropoxy)-3,5-dihydroxy-6-(hydroxymethyl)oxan-4-yl]oxy-6-(hydroxymethyl)-5-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-4-yl]oxy-5-hydroxy-6-(hydroxymethyl)-3-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-4-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]acetamide',
                                      'reason': 'Contains carbonyl group(s)'},
                                  {   'smiles': 'O1[C@@H]([C@H](O)[C@H](O)[C@@H](O)[C@@H]1OC[C@H](OC(=O)C/C=C\\CC)COC(=O)CCCCCCCCCCC)CO[C@H]2O[C@@H]([C@H](O)[C@H](O)[C@H]2O)CO',
                                      'name': 'Capsoside A',
                                      'reason': 'Contains carbonyl group(s)'},
                                  {   'smiles': 'C[C@@H]1O[C@@H](OC[C@@H](O)[C@@H](O)[C@H](O)[C@H](CO)NC(C)=O)[C@@H](O)[C@H](O)[C@@H]1O',
                                      'name': 'alpha-L-Fucp-(1->6)-D-Glc-OH',
                                      'reason': 'Contains carbonyl group(s)'},
                                  {   'smiles': 'P(OC[C@H]1O[C@@H](OC[C@H]2O[C@@H](OC[C@H](OC(=O)CCCCCCC/C=C\\CCCCCC)CO/C=C\\CCCCCCCCCCCC)C(O)C(O)[C@@H]2O)C(O)C(O)[C@@H]1O)(OCCN)(O)=O',
                                      'name': 'PE-GlcGlcDG(P-14:0/16:1(9Z))',
                                      'reason': 'Contains carbonyl group(s)'},
                                  {   'smiles': 'O([C@H]1[C@@H](O)[C@H](O[C@@H](O[C@@H]([C@H](O)[C@@H](O)CO)[C@H](O)CO)[C@@H]1O)CO)[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2NC(=O)C)CO',
                                      'name': 'beta-D-GlcpNAc-(1->3)-beta-D-Galp-(1->4)-D-Glc-OH',
                                      'reason': 'Contains carbonyl group(s)'}],
    'attempt': 0,
    'reasoning': 'Let me create a program to classify alditols. Based on the '
                 "definition and examples, I'll identify key characteristics:\n"
                 '\n'
                 '1. Alditols are acyclic polyols (multiple hydroxyl groups)\n'
                 '2. They have a general formula HOCH2[CH(OH)]nCH2OH\n'
                 '3. They are derived from aldoses by reduction of the '
                 'carbonyl group\n'
                 '4. They typically have multiple chiral centers\n'
                 '5. Key features include:\n'
                 '   - Two terminal CH2OH groups\n'
                 '   - A chain of carbons with OH groups attached\n'
                 '   - No carbonyl groups (C=O)\n'
                 '   - No cyclic structures (for pure alditols)\n'
                 '\n'
                 "Here's the implementation:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 68,
    'num_false_positives': 100,
    'num_true_negatives': 3163,
    'num_false_negatives': 6,
    'num_negatives': None,
    'precision': 0.40476190476190477,
    'recall': 0.918918918918919,
    'f1': 0.5619834710743802,
    'accuracy': 0.9682349415642793,
    'negative_predictive_value': 0.9981066582518144}