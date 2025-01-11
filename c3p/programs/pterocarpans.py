"""
Classifies: CHEBI:26377 pterocarpans
"""
"""
Classifies: pterocarpans
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_pterocarpans(smiles: str):
    """
    Determines if a molecule is a pterocarpan based on its SMILES string.
    Pterocarpans have a characteristic 6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pterocarpan, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic pterocarpan skeleton SMARTS pattern
    # Represents the core 6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene structure
    pterocarpan_core = Chem.MolFromSmarts(
        "[#6]1-[#6]2-[#6]-[#6]3=,:[#6]-[#6](=,:[#6]-[#6]=,:[#6]3-[#8]-[#6]2)-[#8]-[#6]1"
    )
    
    if not mol.HasSubstructMatch(pterocarpan_core):
        return False, "Missing characteristic pterocarpan core structure"

    # Check for the presence of two benzene rings
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    benzene_matches = len(mol.GetSubstructMatches(benzene_pattern))
    if benzene_matches < 2:
        return False, "Missing required benzene rings"

    # Check for the characteristic oxygen bridges (furan and pyran rings)
    oxygen_bridge_pattern = Chem.MolFromSmarts("[#6]~[#8]~[#6]")
    oxygen_bridges = len(mol.GetSubstructMatches(oxygen_bridge_pattern))
    if oxygen_bridges < 2:
        return False, "Missing characteristic oxygen bridges"

    # Check for sp3 carbons at 6a and 11a positions
    sp3_carbon_pattern = Chem.MolFromSmarts("[CX4]")
    sp3_carbons = len(mol.GetSubstructMatches(sp3_carbon_pattern))
    if sp3_carbons < 2:
        return False, "Missing required sp3 carbons at 6a and 11a positions"

    # Verify the molecule is not a coumestan (which would have a C=O group)
    carbonyl_pattern = Chem.MolFromSmarts("[#6]=O")
    if mol.HasSubstructMatch(carbonyl_pattern):
        return False, "Contains carbonyl group (appears to be a coumestan rather than pterocarpan)"

    # Additional check for the complete ring system
    tetracyclic_system = Chem.MolFromSmarts(
        "[#6]1-2-[#6]-[#6]3=,:[#6]-[#6](=,:[#6]-[#6]=,:[#6]3-[#8]-[#6]-2)-[#8]-[#6]-1"
    )
    if not mol.HasSubstructMatch(tetracyclic_system):
        return False, "Missing complete tetracyclic ring system"

    return True, "Contains characteristic 6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene skeleton"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26377',
                          'name': 'pterocarpans',
                          'definition': 'Members of the class of '
                                        'benzofurochromene with a '
                                        '6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene '
                                        'skeleton and its substituted '
                                        'derivatives. They generally bear '
                                        'structural resemblance to '
                                        'isoflavanoids that possess antibiotic '
                                        'activity and are produced by plant '
                                        'tissues in response to infection. '
                                        'They are the 3,4-dihydroderivatives '
                                        'of coumestans.',
                          'parents': ['CHEBI:38834', 'CHEBI:72544'],
                          'xrefs': [   'MetaCyc:PTEROCARPANS',
                                       'Wikipedia:Pterocarpan'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'Cc1cc(C)c2Cc3c(C)cc(C)c(Cc4c(C)cc(C)c(Cc5c(C)cc(C)c(Cc1c2C)c5C)c4C)c3C',
                                     'name': '4,6,10,12,16,18,22,24,25,26,27,28-dodecamethylcalix[4]arene',
                                     'reason': 'Missing characteristic '
                                               'pterocarpan core structure'},
                                 {   'smiles': 'O=C1C(=C2C=C3[C@]([C@@H](C(C)C)[C@@H]([C@H]3O)OC(=O)C)(C)CC[C@]2(C)CC1)COC(=O)C',
                                     'name': 'Dahliane E',
                                     'reason': 'Missing characteristic '
                                               'pterocarpan core structure'},
                                 {   'smiles': 'O[C@@H]([C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)N[C@@H](CC(C)C)C(O)=O)C',
                                     'name': 'Glu-Thr-Leu',
                                     'reason': 'Missing characteristic '
                                               'pterocarpan core structure'},
                                 {   'smiles': 'CCOc1ccc(NC(=O)C(C)O)cc1',
                                     'name': 'p-Lactophenetide',
                                     'reason': 'Missing characteristic '
                                               'pterocarpan core structure'},
                                 {   'smiles': 'O=C1NCC=2C1=C3C(N([C@@H]4O[C@H]([C@@H](O)[C@H]([C@H]4OC)O)C)C5=C3C=CC=C5)=C6NC7=C(C26)C=CC=C7',
                                     'name': "3'-epi-5'-methoxy-K252d",
                                     'reason': 'Missing characteristic '
                                               'pterocarpan core structure'},
                                 {   'smiles': 'O=C(OC1C(O)C(OC(C1O)CO)OCC2OC(OCC(OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)C(O)C(C2O)O)CCCCCCCCCCCCCCC',
                                     'name': '[(2S)-2-hexadecanoyloxy-3-[(2S,3R,4S,5S,6R)-6-[[(2S,3R,4S,5S,6R)-4-hexadecanoyloxy-3,5-dihydroxy-6-(hydroxymethyl)tetrahydropyran-2-yl]oxymethyl]-3,4,5-trihydroxy-tetrahydropyran-2-yl]oxy-propyl] '
                                             '(9E,12E,15E)-octadeca-9,12,15-trienoate',
                                     'reason': 'Missing characteristic '
                                               'pterocarpan core structure'},
                                 {   'smiles': 'O=C1C2=C(O)C=C(OC)C=C2C(=O)C3=C1[C@@H]([C@@H](O)[C@]([C@@H]3O)(O)C)C',
                                     'name': 'Altersolanol G',
                                     'reason': 'Missing characteristic '
                                               'pterocarpan core structure'},
                                 {   'smiles': '[H][C@]1(O[C@](O)(C[C@H](O)[C@H]1NC(=O)CO)C([O-])=O)[C@H](O)[C@H](O)CO',
                                     'name': 'N-glycoloyl-alpha-neuraminate',
                                     'reason': 'Missing characteristic '
                                               'pterocarpan core structure'},
                                 {   'smiles': 'OC(C(O)C/C=C\\C/C=C\\C/C=C\\CC)C/C=C\\C/C=C\\CCC(O)=O',
                                     'name': '10,11-DiHDPE',
                                     'reason': 'Missing characteristic '
                                               'pterocarpan core structure'},
                                 {   'smiles': '[Na+].[H][C@]12SCC(C)=C(N1C(=O)[C@H]2NC(=O)[C@H](N)c1ccccc1)C([O-])=O',
                                     'name': 'cephalexin sodium',
                                     'reason': 'Missing characteristic '
                                               'pterocarpan core structure'}],
    'sample_false_negatives': [   {   'smiles': '[H][C@@]12COc3cc(O)ccc3[C@]1([H])Oc1c(CCC(C)(C)O)c(O)ccc21',
                                      'name': 'phaseollidin hydrate',
                                      'reason': 'Missing characteristic '
                                                'pterocarpan core structure'},
                                  {   'smiles': 'COc1ccc2[C@@H]3COc4cc(O[C@@H]5O[C@H](COC(=O)CC(O)=O)[C@@H](O)[C@H](O)[C@H]5O)ccc4[C@@H]3Oc2c1',
                                      'name': 'Medicarpin '
                                              "3-O-glucoside-6'-malonate",
                                      'reason': 'Missing characteristic '
                                                'pterocarpan core structure'},
                                  {   'smiles': '[H][C@]12Oc3cc4OCOc4cc3[C@@]1(O)COc1cc(O)ccc21',
                                      'name': '(+)-6a-hydroxymaackiain',
                                      'reason': 'Missing characteristic '
                                                'pterocarpan core structure'},
                                  {   'smiles': 'O1C2C(C3=C1C(C[C@H](O)C(C)=C)=C(O)C=C3)COC4=C2C=CC(O)=C4',
                                      'name': 'Dolichin B',
                                      'reason': 'Missing characteristic '
                                                'pterocarpan core structure'},
                                  {   'smiles': 'O1C2C(C3=C1C=C(OC)C=C3)COC4=C2C=CC(O)=C4OC',
                                      'name': '3-Hydroxy-4,9-dimethoxypterocarpan',
                                      'reason': 'Missing characteristic '
                                                'pterocarpan core structure'},
                                  {   'smiles': 'O1[C@@]2([C@](C3=C1C=C4OC(C=CC4=C3)(C)C)(COC5=C2C(OC)=C(C(O)=C5)CC=C(C)C)[H])[H]',
                                      'name': 'Kanzonol F',
                                      'reason': 'Missing characteristic '
                                                'pterocarpan core structure'},
                                  {   'smiles': 'O[C@@]12COc3c4OCOc4ccc3[C@@H]1Oc1cc3OCOc3cc21',
                                      'name': '(-)-Acanthocarpan',
                                      'reason': 'Missing characteristic '
                                                'pterocarpan core structure'},
                                  {   'smiles': 'COc1c(CC=C(C)C)c(O)cc2OC[C@@H]3[C@@H](Oc4cc(O)ccc34)c12',
                                      'name': 'edudiol',
                                      'reason': 'Missing characteristic '
                                                'pterocarpan core structure'},
                                  {   'smiles': '[H][C@@]12COc3cc(OC)ccc3[C@]1([H])Oc1c(CC=C(C)C)c(OC)c(O)cc21',
                                      'name': 'lespeflorin G8',
                                      'reason': 'Missing characteristic '
                                                'pterocarpan core structure'},
                                  {   'smiles': '[H][C@]12COc3ccccc3[C@@]1([H])Oc1ccccc21',
                                      'name': '(6aS,11aS)-pterocarpan',
                                      'reason': 'Missing characteristic '
                                                'pterocarpan core structure'}],
    'attempt': 0,
    'reasoning': "I'll help create a program to classify pterocarpans. Let's "
                 'break this down:\n'
                 '\n'
                 'Key characteristics of pterocarpans:\n'
                 '1. They have a benzofuro[3,2-c]chromene skeleton\n'
                 '2. They are 6a,11a-dihydro derivatives\n'
                 '3. They contain a specific tetracyclic system with two '
                 'benzene rings connected through a dihydrofuran and a '
                 'dihydropyran ring\n'
                 '4. They are 3,4-dihydro derivatives of coumestans\n'
                 '\n'
                 "Here's the program:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 142212,
    'num_false_negatives': 88,
    'num_negatives': None,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': 0.9993815881939564,
    'negative_predictive_value': 0.9993815881939564}