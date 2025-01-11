"""
Classifies: CHEBI:2440 acrovestone
"""
"""
Classifies: acrovestone compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_acrovestone(smiles: str):
    """
    Determines if a molecule is an acrovestone compound based on its SMILES string.
    Acrovestone compounds are typically isoflavone derivatives with glycosidic attachments.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acrovestone compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for chromone core (benzopyran with ketone)
    chromone_pattern = Chem.MolFromSmarts("O=C1C=COc2ccccc12")
    if not mol.HasSubstructMatch(chromone_pattern):
        return False, "No chromone core structure found"

    # Check for presence of sugar moiety patterns
    sugar_pattern = Chem.MolFromSmarts("[OX2][CH]1[CH][CH][CH]([CH][CH]1)(O)[CH2][OX2]")
    if not mol.HasSubstructMatch(sugar_pattern):
        sugar_pattern2 = Chem.MolFromSmarts("[OX2][CH]1[CH][CH][CH]([CH]1)(O)[CH2][OX2]")
        if not mol.HasSubstructMatch(sugar_pattern2):
            return False, "No glycosidic moiety found"

    # Count hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if hydroxyl_matches < 2:
        return False, "Insufficient hydroxyl groups"

    # Check for aromatic rings
    aromatic_rings = len(rdMolDescriptors.CalcAromaticRings(mol))
    if aromatic_rings < 2:
        return False, "Insufficient aromatic rings"

    # Calculate molecular weight - should be substantial due to sugar moieties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for acrovestone compound"

    # Count oxygen atoms - should have multiple due to glycosidic bonds and hydroxyl groups
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 5:
        return False, "Insufficient oxygen atoms"

    # Look for potential methoxy groups (common in these compounds)
    methoxy_pattern = Chem.MolFromSmarts("COc")
    methoxy_matches = len(mol.GetSubstructMatches(methoxy_pattern))

    # Build reason string
    features = []
    features.append("Contains chromone core")
    features.append(f"Has {hydroxyl_matches} hydroxyl groups")
    features.append(f"Contains {aromatic_rings} aromatic rings")
    features.append("Contains glycosidic moiety")
    if methoxy_matches > 0:
        features.append(f"Has {methoxy_matches} methoxy groups")
    
    return True, "; ".join(features)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:2440',
                          'name': 'acrovestone',
                          'definition': 'A polyphenol that is isolated from '
                                        'Acronychia pedunculata and exhibits '
                                        'moderate antioxidant and '
                                        'antityrosinase activities.',
                          'parents': [   'CHEBI:22187',
                                         'CHEBI:26195',
                                         'CHEBI:35618',
                                         'CHEBI:78840'],
                          'xrefs': [   'CAS:24177-16-0',
                                       'KEGG:C09916',
                                       'KNApSAcK:C00002978',
                                       'PMID:12880321',
                                       'PMID:22708987',
                                       'PMID:2506892',
                                       'PMID:2614422',
                                       'Reaxys:2611012'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'C[C@H]1CN(C(=O)C2=C(C3=CC=CC=C3CO[C@@H]1CN(C)C(=O)NC4=CC=C(C=C4)OC)C5=CC=CC=C5N2C)[C@H](C)CO',
                                     'name': 'LSM-30659',
                                     'reason': 'No chromone core structure '
                                               'found'},
                                 {   'smiles': 'c1c[nH+]c[nH]1',
                                     'name': 'imidazolium cation',
                                     'reason': 'No chromone core structure '
                                               'found'},
                                 {   'smiles': 'ClC=1C(O)=C2O[C@@]([C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC)[C@@H](O)C(=C2)C1)C(C)=C)C(=O)N3[C@@H](C=CC3)C(=O)N/C(=C(/CC)\\C)/C(=O)N/C(=C/C(O)=O)/C(O)=O)(CC)C',
                                     'name': 'Phomopsin A',
                                     'reason': 'No chromone core structure '
                                               'found'},
                                 {   'smiles': 'ClC(Cl)[C@H](O)CC=1OC(=O)C=2C(O)=CC(=CC2C1)O',
                                     'name': 'Desmethyldichlorodiaportin',
                                     'reason': 'No chromone core structure '
                                               'found'},
                                 {   'smiles': 'ClC=1C(=C(O)C2=C(C1O)C(=O)C=CC2=O)CC=C(C)C',
                                     'name': 'Chlorosesamone',
                                     'reason': 'No chromone core structure '
                                               'found'},
                                 {   'smiles': 'CCCCC[C@H]1O[C@H]1C\\C=C/CCCCCCCC(O)=O',
                                     'name': '(+)-vernolic acid',
                                     'reason': 'No chromone core structure '
                                               'found'},
                                 {   'smiles': 'O(C1C(O)C(OC(OC=2C3=C(C=CC2)C=C(C(=C3O)C(=O)C)C)C1O)CO)C4OC(C(O)C(O)C4O)CO',
                                     'name': 'Orientaloside',
                                     'reason': 'No chromone core structure '
                                               'found'},
                                 {   'smiles': 'O(C1=C(C=CC=C1C(O)=O)C)C(=O)C',
                                     'name': 'CRESOPYRINE',
                                     'reason': 'No chromone core structure '
                                               'found'},
                                 {   'smiles': 'C[C@H]1CN(S(=O)(=O)C2=C(C=C(C=C2)C#CC3=CC=NC=C3)O[C@@H]1CN(C)C(=O)NC4=CC5=C(C=C4)OCO5)[C@@H](C)CO',
                                     'name': '3-(1,3-benzodioxol-5-yl)-1-[[(4S,5S)-2-[(2S)-1-hydroxypropan-2-yl]-4-methyl-1,1-dioxo-8-(2-pyridin-4-ylethynyl)-4,5-dihydro-3H-6,1$l^{6},2-benzoxathiazocin-5-yl]methyl]-1-methylurea',
                                     'reason': 'No chromone core structure '
                                               'found'},
                                 {   'smiles': 'O1C2C3C(CCC3=C)C(CCC2C(C1=O)=C)=C',
                                     'name': '3,6,9-Trimethylidene-3a,4,5,6a,7,8,9a,9b-octahydroazuleno[4,5-b]furan-2-one',
                                     'reason': 'No chromone core structure '
                                               'found'}],
    'sample_false_negatives': [   {   'smiles': 'O1C([C@H](O)C(O)C(O)[C@H]1OC2=CC=C(C=C2)C=3C(=O)C=4C(OC3)=CC(O)=CC4O)CO',
                                      'name': "Genistein 4'-O-glucoside",
                                      'reason': 'No chromone core structure '
                                                'found'},
                                  {   'smiles': 'S(OC1=CC=C(C=2C(=O)C3=C(OC2)C=C(O[C@@H]4O[C@@H]([C@@H](O)[C@H](O)[C@H]4O)C(O)=O)C=C3)C=C1)(O)(=O)=O',
                                      'name': 'daidzein '
                                              "7-glucuronide-4'-sulfate",
                                      'reason': 'No chromone core structure '
                                                'found'},
                                  {   'smiles': 'O1[C@@H]([C@@H](O)[C@H](O)[C@@H](O)[C@@H]1OC2=CC=3OC=C(C(=O)C3C=C2)C4=CC=5OCOC5C=C4)CO',
                                      'name': 'Rothindin',
                                      'reason': 'No chromone core structure '
                                                'found'},
                                  {   'smiles': 'O1[C@@H]([C@@H](O)[C@H](O)[C@@H](O)[C@@H]1OC2=CC=3OC=C(C4=CC=C(O[C@@H]5O[C@@H]([C@@H](O)[C@H](O)[C@H]5O)C(O)=O)C=C4)C(=O)C3C=C2)CO',
                                      'name': "Daidzin 4'-O-glucuronide",
                                      'reason': 'No chromone core structure '
                                                'found'},
                                  {   'smiles': 'O1C(C(O)C(O)C(O)C1OC2=CC=3OC=C(C(=O)C3C=C2)C4=CC=C(OC)C=C4)COC(=O)CC(O)=O',
                                      'name': 'Formononetin '
                                              "7-(6''-malonylglucoside)",
                                      'reason': 'No chromone core structure '
                                                'found'},
                                  {   'smiles': 'O1C([C@@H](O)[C@H](O)C(O)[C@@H]1OC2=CC=3OCC(CC3C=C2)C4=C(O)C(OC)=C(OC)C=C4)CO',
                                      'name': 'Isomucronulator 7-O-glucoside',
                                      'reason': 'No chromone core structure '
                                                'found'},
                                  {   'smiles': 'O1C(C(O)C(O)C(O)C1OC2=CC=C(C=C2)C=3C(=O)C=4C(OC3)=CC=5OCOC5C4O)CO',
                                      'name': "Irilone 4'-glucoside",
                                      'reason': 'No chromone core structure '
                                                'found'},
                                  {   'smiles': 'O(C1=CC2=C(C(C(=CO2)C3=CC(=C(C=C3OC)OC)OC)=O)C=C1)[C@H]4[C@H](C([C@@H](C(O4)CO)O)O)O',
                                      'name': "7-hydroxy-2',4',5'-trimethoxyisoflavone "
                                              '7-O-glucoside',
                                      'reason': 'No chromone core structure '
                                                'found'},
                                  {   'smiles': 'O1C(C(O)C(O)C(O)C1OC2=CC=3OC=C(C(=O)C3C=C2)C4=CC=C(O)C=C4)COC(=O)C',
                                      'name': '[6-[3-(4-Hydroxyphenyl)-4-oxidanylidene-chromen-7-yl]oxy-3,4,5-tris(oxidanyl)oxan-2-yl]methyl '
                                              'ethanoate',
                                      'reason': 'No chromone core structure '
                                                'found'},
                                  {   'smiles': 'O1C(C(O)C(O)C(O)C1OC=2C=C3OC=C(C(=O)C3=C(O)C2)C=4C(O)=CC=CC4)CO',
                                      'name': 'Isogenistein 7-glucoside',
                                      'reason': 'No chromone core structure '
                                                'found'}],
    'attempt': 0,
    'reasoning': 'Let me create a program to classify acrovestone compounds. '
                 'From analyzing the examples, I notice that these compounds '
                 'are typically isoflavone derivatives with various glycoside '
                 'attachments. The key structural features appear to be:\n'
                 '\n'
                 '1. A chromone core structure (benzopyran with a ketone)\n'
                 '2. Often have glycosidic linkages (sugar moieties)\n'
                 '3. Multiple phenolic/hydroxyl groups\n'
                 '4. Various substitution patterns including methoxy groups\n'
                 '5. Sometimes contain additional ring systems\n'
                 '\n'
                 "Here's the program to identify acrovestone compounds:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 142215,
    'num_false_negatives': 85,
    'num_negatives': None,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': 0.9994026704146171,
    'negative_predictive_value': 0.9994026704146171}