"""
Classifies: CHEBI:83411 beta-D-glucosiduronate
"""
"""
Classifies: beta-D-glucosiduronate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_D_glucosiduronate(smiles: str):
    """
    Determines if a molecule contains a beta-D-glucosiduronate group.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule contains beta-D-glucosiduronate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for beta-D-glucuronic acid
    # Pattern matches:
    # - O-glycosidic bond
    # - Beta configuration at C1
    # - Correct stereochemistry at C2,C3,C4
    # - Carboxylate group at C5
    glucuronate_pattern = Chem.MolFromSmarts(
        "[OX2][C@@H]1[C@@H]([C@H]([C@@H]([C@H](O1)C(=O)[O-])O)O)O"
    )
    
    if not mol.HasSubstructMatch(glucuronate_pattern):
        return False, "No beta-D-glucuronate group found with correct stereochemistry"
    
    # Count matches to ensure we have at least one
    matches = mol.GetSubstructMatches(glucuronate_pattern)
    if len(matches) < 1:
        return False, "No beta-D-glucuronate group found"
        
    # Verify presence of carboxylate
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)[O-]")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylate group found"
    
    # Check for required hydroxyls
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyls = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if hydroxyls < 3:
        return False, "Insufficient hydroxyl groups"
        
    # Verify glycosidic bond
    glycosidic_pattern = Chem.MolFromSmarts("[OX2][CH]1O[CH][CH][CH][CH]1")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No glycosidic bond found"

    return True, "Contains beta-D-glucuronate group with correct stereochemistry and carboxylate"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:83411',
                          'name': 'beta-D-glucosiduronate',
                          'definition': 'A carbohydrate acid derivative anion '
                                        'obtained by deprotonation of the '
                                        'carboxy group of any '
                                        'beta-D-glucosiduronic acid; major '
                                        'species at pH 7.3.',
                          'parents': ['CHEBI:132367', 'CHEBI:63551'],
                          'xrefs': ['MetaCyc:Beta-D-Glucuronides'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'O1[C@H](CC2=C(C=3C=4C(C(O)=C5C3C[C@H](C)OC5)=C(OC)C=C(C4)OC)C6=CC(OC)=CC(=C6C(=C2C1)O)OC)C',
                                     'name': '(3S)-5-[(3S)-10-hydroxy-7,9-dimethoxy-3-methyl-3,4-dihydro-1H-benzo[g]isochromen-5-yl]-7,9-dimethoxy-3-methyl-3,4-dihydro-1H-benzo[g]isochromen-10-ol',
                                     'reason': 'No beta-D-glucuronate group '
                                               'found with correct '
                                               'stereochemistry'},
                                 {   'smiles': 'C[C@H]1CN([C@@H](COC2=C(C=CC(=C2)NC(=O)C)C(=O)N(C[C@@H]1OC)C)C)C(=O)C3CCC3',
                                     'name': 'N-[(5R,6S,9R)-8-[cyclobutyl(oxo)methyl]-5-methoxy-3,6,9-trimethyl-2-oxo-11-oxa-3,8-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]acetamide',
                                     'reason': 'No beta-D-glucuronate group '
                                               'found with correct '
                                               'stereochemistry'},
                                 {   'smiles': '[C@@]12(CCCC[C@@]1(C=C[C@@H]([C@@H]2CC[C@H](C[C@H](CC(=O)[O-])O)O)C)[H])[H]',
                                     'name': '4a,5-dihydro-ML-236C carboxylate',
                                     'reason': 'No beta-D-glucuronate group '
                                               'found with correct '
                                               'stereochemistry'},
                                 {   'smiles': 'S(=O)(CC1OC(N2C3=NC=NC(N)=C3N=C2)C(O)C1O)C',
                                     'name': "(S)-5'-Deoxy-5'-(methylsulfinyl)adenosine",
                                     'reason': 'No beta-D-glucuronate group '
                                               'found with correct '
                                               'stereochemistry'},
                                 {   'smiles': 'OC=1C(=C(C=2C=3C(NC2)=CC=CC3)C(=O)C(=O)C1C=4C=5C(NC4)=CC=CC5)CC=C(C)C',
                                     'name': 'Ochrindole D',
                                     'reason': 'No beta-D-glucuronate group '
                                               'found with correct '
                                               'stereochemistry'},
                                 {   'smiles': 'O(C1C(O)C(OC1OC2=C(OC=3C(C2=O)=C(O)C=C(OC4OC(C(O)C(O)C4O)C)C3)C5=CC=C(O)C=C5)CO)C6OCC(O)(C6O)CO',
                                     'name': 'Kaempferol '
                                             '3-apiosyl-(1->2)-alpha-L-arabinofuranoside-7-rhamnoside',
                                     'reason': 'No beta-D-glucuronate group '
                                               'found with correct '
                                               'stereochemistry'},
                                 {   'smiles': '[O-]C(=O)OON=O',
                                     'name': 'nitrosoperoxycarbonate(1-)',
                                     'reason': 'No beta-D-glucuronate group '
                                               'found with correct '
                                               'stereochemistry'},
                                 {   'smiles': 'P(OCC(O)COC(=O)CCCCCCCCCCCCCCCCC)(O)(O)=O',
                                     'name': 'LysoPA(18:0/0:0)',
                                     'reason': 'No beta-D-glucuronate group '
                                               'found with correct '
                                               'stereochemistry'},
                                 {   'smiles': 'O=C1C=2C(=O)OC34C2OC(C(O)C3C=CC(C(CC1C)C)=O)C(C)C4',
                                     'name': 'Atrop-Abybetaomicin C',
                                     'reason': 'No beta-D-glucuronate group '
                                               'found with correct '
                                               'stereochemistry'},
                                 {   'smiles': 'P(OC[C@H](OC(=O)CCCCCCCCCCCCCCCCCC)CO/C=C\\CCCCCCCCCCCCCCCC)(OC[C@@H](O)CO)(O)=O',
                                     'name': 'PG(P-18:0/19:0)',
                                     'reason': 'No beta-D-glucuronate group '
                                               'found with correct '
                                               'stereochemistry'}],
    'sample_false_negatives': [],
    'attempt': 0,
    'reasoning': 'Let me help create a program to classify '
                 'beta-D-glucosiduronates. Based on the examples, these '
                 'compounds have some key characteristics:\n'
                 '\n'
                 '1. They contain a beta-D-glucuronic acid moiety attached via '
                 'a glycosidic bond\n'
                 '2. The glucuronic acid has a carboxylate group ([O-])\n'
                 '3. The glucuronic acid follows a specific stereochemistry '
                 'pattern\n'
                 '4. The glycosidic bond is in beta configuration\n'
                 '\n'
                 "Here's the program:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 44,
    'num_false_positives': 100,
    'num_true_negatives': 132650,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.3055555555555556,
    'recall': 1.0,
    'f1': 0.46808510638297873,
    'accuracy': 0.9992469539286414,
    'negative_predictive_value': 1.0}