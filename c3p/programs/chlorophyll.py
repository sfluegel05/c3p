"""
Classifies: CHEBI:28966 chlorophyll
"""
"""
Classifies: chlorophyll molecules
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_chlorophyll(smiles: str):
    """
    Determines if a molecule is a chlorophyll based on its SMILES string.
    Chlorophylls are magnesium porphyrins with a fifth ring and various side chains.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chlorophyll, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for magnesium atom (including different charge states)
    mg_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'Mg']
    if not mg_atoms:
        return False, "No magnesium atom found"
    if len(mg_atoms) > 1:
        return False, "More than one magnesium atom found"

    # Check for 4 nitrogens in the molecule
    n_atoms = len([atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'N'])
    if n_atoms != 4:
        return False, f"Must have exactly 4 nitrogens, found {n_atoms}"

    # Simplified porphyrin-like core pattern
    porphyrin_pattern = Chem.MolFromSmarts("[#7]~[#6]~[#6]~[#6]~[#7]~[#6]~[#6]~[#6]~[#7]~[#6]~[#6]~[#6]~[#7]")
    if not mol.HasSubstructMatch(porphyrin_pattern):
        return False, "No porphyrin-like core found"

    # Check for fifth ring (cyclopentanone or similar)
    fifth_ring_pattern = Chem.MolFromSmarts("[#6]1~[#6]~[#6]~[#6]~[#6]1")
    if not mol.HasSubstructMatch(fifth_ring_pattern):
        return False, "No characteristic fifth ring found"

    # Count rings
    rings = mol.GetRingInfo()
    if rings.NumRings() < 5:
        return False, "Too few rings for chlorophyll structure"

    # Check for characteristic substituents
    substituents = []
    
    # Vinyl group (-CH=CH2)
    vinyl_pattern = Chem.MolFromSmarts("C=C")
    if mol.HasSubstructMatch(vinyl_pattern):
        substituents.append("vinyl")
        
    # Ethyl group
    ethyl_pattern = Chem.MolFromSmarts("CC")
    if mol.HasSubstructMatch(ethyl_pattern):
        substituents.append("ethyl")
    
    # Carboxyl/ester groups
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if mol.HasSubstructMatch(carboxyl_pattern):
        substituents.append("carboxyl/ester")

    # Long chain (phytol or similar)
    long_chain = Chem.MolFromSmarts("CCCC")
    if mol.HasSubstructMatch(long_chain):
        substituents.append("alkyl chain")

    if not substituents:
        return False, "Missing characteristic substituents"

    # Count carbons and check molecular complexity
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, "Too few carbons for chlorophyll structure"

    # Additional check for conjugated system
    conjugated_pattern = Chem.MolFromSmarts("C=C-C=C")
    if not mol.HasSubstructMatch(conjugated_pattern):
        return False, "Missing conjugated system"

    return True, f"Contains magnesium-coordinated porphyrin core with fifth ring and {', '.join(substituents)}"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:28966',
                          'name': 'chlorophyll',
                          'definition': 'A family of magnesium porphyrins, '
                                        'defined by the presence of a fifth '
                                        'ring beyond the four pyrrole-like '
                                        'rings. The rings can have various '
                                        'side chains which usually include a '
                                        'long phytol chain.',
                          'parents': ['CHEBI:25111'],
                          'xrefs': [   'CAS:1406-65-1',
                                       'COMe:MOL000012',
                                       'KEGG:C01793',
                                       'PMID:29286160'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Error: Python argument types in\n'
               '    Mol.HasSubstructMatch(Mol, NoneType)\n'
               'did not match C++ signature:\n'
               '    HasSubstructMatch(RDKit::ROMol self, RDKit::MolBundle '
               'query, RDKit::SubstructMatchParameters params=True)\n'
               '    HasSubstructMatch(RDKit::ROMol self, RDKit::ROMol query, '
               'RDKit::SubstructMatchParameters params)\n'
               '    HasSubstructMatch(RDKit::ROMol self, RDKit::MolBundle '
               'query, bool recursionPossible=True, bool useChirality=False, '
               'bool useQueryQueryMatches=False)\n'
               '    HasSubstructMatch(RDKit::ROMol self, RDKit::ROMol query, '
               'bool recursionPossible=True, bool useChirality=False, bool '
               'useQueryQueryMatches=False)\n'
               'Attempt failed: F1 score of 0 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: NONE\n'
               'False positives: NONE\n'
               'False negatives: NONE\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'Cc1cc(C)c2Cc3c(C)cc(C)c(Cc4c(C)cc(C)c(Cc5c(C)cc(C)c(Cc1c2C)c5C)c4C)c3C',
                                     'name': '4,6,10,12,16,18,22,24,25,26,27,28-dodecamethylcalix[4]arene',
                                     'reason': 'No magnesium atom found'},
                                 {   'smiles': 'O=C1C(=C2C=C3[C@]([C@@H](C(C)C)[C@@H]([C@H]3O)OC(=O)C)(C)CC[C@]2(C)CC1)COC(=O)C',
                                     'name': 'Dahliane E',
                                     'reason': 'No magnesium atom found'},
                                 {   'smiles': 'O[C@@H]([C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)N[C@@H](CC(C)C)C(O)=O)C',
                                     'name': 'Glu-Thr-Leu',
                                     'reason': 'No magnesium atom found'},
                                 {   'smiles': 'CCOc1ccc(NC(=O)C(C)O)cc1',
                                     'name': 'p-Lactophenetide',
                                     'reason': 'No magnesium atom found'},
                                 {   'smiles': 'O=C1NCC=2C1=C3C(N([C@@H]4O[C@H]([C@@H](O)[C@H]([C@H]4OC)O)C)C5=C3C=CC=C5)=C6NC7=C(C26)C=CC=C7',
                                     'name': "3'-epi-5'-methoxy-K252d",
                                     'reason': 'No magnesium atom found'},
                                 {   'smiles': 'O=C(OC1C(O)C(OC(C1O)CO)OCC2OC(OCC(OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)C(O)C(C2O)O)CCCCCCCCCCCCCCC',
                                     'name': '[(2S)-2-hexadecanoyloxy-3-[(2S,3R,4S,5S,6R)-6-[[(2S,3R,4S,5S,6R)-4-hexadecanoyloxy-3,5-dihydroxy-6-(hydroxymethyl)tetrahydropyran-2-yl]oxymethyl]-3,4,5-trihydroxy-tetrahydropyran-2-yl]oxy-propyl] '
                                             '(9E,12E,15E)-octadeca-9,12,15-trienoate',
                                     'reason': 'No magnesium atom found'},
                                 {   'smiles': 'O=C1C2=C(O)C=C(OC)C=C2C(=O)C3=C1[C@@H]([C@@H](O)[C@]([C@@H]3O)(O)C)C',
                                     'name': 'Altersolanol G',
                                     'reason': 'No magnesium atom found'},
                                 {   'smiles': '[H][C@]1(O[C@](O)(C[C@H](O)[C@H]1NC(=O)CO)C([O-])=O)[C@H](O)[C@H](O)CO',
                                     'name': 'N-glycoloyl-alpha-neuraminate',
                                     'reason': 'No magnesium atom found'},
                                 {   'smiles': 'OC(C(O)C/C=C\\C/C=C\\C/C=C\\CC)C/C=C\\C/C=C\\CCC(O)=O',
                                     'name': '10,11-DiHDPE',
                                     'reason': 'No magnesium atom found'},
                                 {   'smiles': '[Na+].[H][C@]12SCC(C)=C(N1C(=O)[C@H]2NC(=O)[C@H](N)c1ccccc1)C([O-])=O',
                                     'name': 'cephalexin sodium',
                                     'reason': 'No magnesium atom found'}],
    'sample_false_negatives': [   {   'smiles': 'C=12N3C(=CC4=[N+]5C(=CC=6N7C=8C(=C9[N+](=C(C1)[C@H]([C@@H]9CCC(O)=O)C)[Mg-2]735)[C@H](C(C8C6C)=O)C(=O)OC)[C@@H]([C@H]4C)CC)C(=C2C)C(=O)C',
                                      'name': 'bacteriochlorophyllide a',
                                      'reason': 'Missing conjugated system'},
                                  {   'smiles': 'COC(=O)[C@H]1C(=O)c2c(C)c3=CC4=N\\C(=C/c5c(C(C)=O)c(C)c6\\C=C7/N=C([C@@H](CCC(=O)OC\\C=C(/C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)[C@@H]7C)C1=c2n3[Mg]n56)[C@H](C)C/4=C/C',
                                      'name': '(7R,8Z)-bacteriochlorophyll b',
                                      'reason': 'Missing conjugated system'},
                                  {   'smiles': '[H]C(C)=C1C(C)C2=CC3=C(C(C)=O)C(C)=C4C=C5[C@@H](C)[C@H](CCC(=O)OC\\C=C(/C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)C6=[N+]5[Mg--]5(N34)N3C(=CC1=[N+]25)C(C)=C1C(=O)[C@H](C(=O)OC)C6=C31',
                                      'name': 'bacteriochlorophyll b',
                                      'reason': 'Missing conjugated system'},
                                  {   'smiles': 'CC[C@@H]1[C@@H](C)C2=Cc3c(C=C)c(C)c4C=C5[C@@H](C)[C@H](CCC(O)=O)C6=[N+]5[Mg--]5(n34)n3c(=CC1=[N+]25)c(C)c1C(=O)[C@H](C(=O)OC)C6=c31',
                                      'name': '3-vinylbacteriochlorophyllide a',
                                      'reason': 'Missing conjugated system'},
                                  {   'smiles': 'CC[C@@H]1[C@@H](C)C2=Cc3c(C(C)=O)c(C)c4C=C5[C@@H](C)[C@H](CCC(=O)OC\\C=C(/C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)C6=[N+]5[Mg--]5(n34)n3c(=CC1=[N+]25)c(C)c1C(=O)[C@H](C(=O)OC)C6=c31',
                                      'name': 'bacteriochlorophyll a',
                                      'reason': 'Missing conjugated system'}],
    'attempt': 3,
    'reasoning': "Let's analyze what went wrong and improve the program:\n"
                 '\n'
                 '1. The main error occurred in the SMARTS pattern matching, '
                 'specifically with `HasSubstructMatch`. The error suggests '
                 "that one of the SMARTS patterns wasn't properly created "
                 '(returned None).\n'
                 '\n'
                 '2. The SMARTS patterns were too complex and included '
                 'newlines and spaces, which can cause parsing issues. We '
                 'should simplify them and remove whitespace.\n'
                 '\n'
                 '3. The porphyrin core pattern might have been too '
                 'restrictive. Chlorophylls can have variations in their core '
                 'structure.\n'
                 '\n'
                 '4. We need to better handle the Mg coordination, as these '
                 'molecules often use special notation like [Mg--] or [Mg-2].\n'
                 '\n'
                 "Here's an improved version:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 22,
    'num_false_positives': 21,
    'num_true_negatives': 142252,
    'num_false_negatives': 5,
    'num_negatives': None,
    'precision': 0.5116279069767442,
    'recall': 0.8148148148148148,
    'f1': 0.6285714285714286,
    'accuracy': 0.9998172874209417,
    'negative_predictive_value': 0.9999648523447001}