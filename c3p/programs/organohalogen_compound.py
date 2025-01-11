"""
Classifies: CHEBI:17792 organohalogen compound
"""
"""
Classifies: organohalogen compound
A compound containing at least one carbon-halogen bond (where X is a halogen atom)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_organohalogen_compound(smiles: str):
    """
    Determines if a molecule is an organohalogen compound based on its SMILES string.
    An organohalogen compound contains at least one carbon-halogen bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organohalogen compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define halogen atomic numbers
    halogens = {9: 'F', 17: 'Cl', 35: 'Br', 53: 'I'}
    
    # Look for halogen atoms in the molecule
    halogen_atoms = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in halogens:
            halogen_atoms.append(atom)
    
    if not halogen_atoms:
        return False, "No halogen atoms found"
    
    # Check if any halogen is bonded to carbon
    c_x_bonds = []  # List to store carbon-halogen bonds
    for halogen in halogen_atoms:
        # Get neighboring atoms
        neighbors = halogen.GetNeighbors()
        for neighbor in neighbors:
            if neighbor.GetAtomicNum() == 6:  # Carbon atomic number
                halogen_symbol = halogens[halogen.GetAtomicNum()]
                c_x_bonds.append(f"C-{halogen_symbol}")
                
    if not c_x_bonds:
        return False, "No carbon-halogen bonds found"
    
    # Create detailed message about the types of C-X bonds found
    unique_bonds = set(c_x_bonds)
    bond_counts = {bond: c_x_bonds.count(bond) for bond in unique_bonds}
    bond_description = ", ".join(f"{count} {bond}" for bond, count in bond_counts.items())
    
    return True, f"Contains carbon-halogen bonds: {bond_description}"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:17792',
                          'name': 'organohalogen compound',
                          'definition': 'A compound containing at least one '
                                        'carbon-halogen bond (where X is a '
                                        'halogen atom).',
                          'parents': ['CHEBI:33285', 'CHEBI:37578'],
                          'xrefs': [   'KEGG:C01322',
                                       'MetaCyc:Organohalogen-Compounds'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'C[C@@H](C1=CC=CC=C1)NC(=O)C[C@H]2CC[C@@H]([C@@H](O2)CO)NC(=O)CN3CCOCC3',
                                     'name': '2-[(2R,5S,6R)-6-(hydroxymethyl)-5-[[2-(4-morpholinyl)-1-oxoethyl]amino]-2-oxanyl]-N-[(1S)-1-phenylethyl]acetamide',
                                     'reason': 'No halogen atoms found'},
                                 {   'smiles': 'C[N+](C)(C)[C@@H](Cc1c[nH]c(n1)S(=O)C[C@H](NC(=O)CC[C@H]([NH3+])C([O-])=O)C([O-])=O)C([O-])=O',
                                     'name': 'N(alpha)-(L-gamma-glutamyl)-hercynyl-L-cysteine '
                                             'sulfoxide(1-)',
                                     'reason': 'No halogen atoms found'},
                                 {   'smiles': 'O(C1=CC=2[C@]3([C@](N(CC3)C)(N(C2C=C1)C)[H])C)C(=O)N4CCC=5C(C4)=CC=CC5',
                                     'name': 'quilostigmine',
                                     'reason': 'No halogen atoms found'},
                                 {   'smiles': 'O[C@@H]1[C@]23[C@@]4(N(C[C@@]([C@]2(C[C@@]4([C@]56[C@]3(CC(=O)[C@](C5)(C([C@H]6O)=C)[H])[H])[H])[H])(CC1)C)CC)[H]',
                                     'name': 'Bullatine G',
                                     'reason': 'No halogen atoms found'},
                                 {   'smiles': 'NC1=NC=NC2=C1N=CN2[C@@H]3O[C@H](COP(=O)(O)O)[C@@H](OC(=O)[C@@H](N)CCC(O)=O)[C@H]3O',
                                     'name': "3'-L-glutamyl-AMP",
                                     'reason': 'No halogen atoms found'},
                                 {   'smiles': 'O1C2(C(C3(C(C4(C(CC3OC(=O)C)C(OC(=O)CC4)(C)C)C)CC2)C)CC15C6N(C=7C5=CC=CC7)C(=O)C(N6)C)C',
                                     'name': 'Teraspiridole C_130091',
                                     'reason': 'No halogen atoms found'},
                                 {   'smiles': 'O(C1[C@@H](OC(=O)C)C(O[C@@H](OC2=C(OC3=C(C2=O)C(O)=CC(O[C@@H]4OC([C@@H](O)[C@H](O)C4O)CO)=C3CC=C(C)C)C5=CC=C(OC)C=C5)[C@H]1O)C)[C@@H]6OC[C@@H](O)[C@H](OC(=O)C)C6O',
                                     'name': 'Sempervirenoside A',
                                     'reason': 'No halogen atoms found'},
                                 {   'smiles': 'COC[C@]1(C(=O)C2CCN1CC2)CO',
                                     'name': '(2S)-2-(hydroxymethyl)-2-(methoxymethyl)-1-azabicyclo[2.2.2]octan-3-one',
                                     'reason': 'No halogen atoms found'},
                                 {   'smiles': 'O[C@H]1CC=2C(N(C=3C1=CC=CC3)C(=O)N)=CC=CC2',
                                     'name': '(S)-MHD',
                                     'reason': 'No halogen atoms found'},
                                 {   'smiles': 'S(OC=1C(O)=C(\\C=C\\C2=CC(O)=C(CC=C(C)C)C(O)=C2)C=CC1O)(O)(=O)=O',
                                     'name': '3-{(e)-2-[3,5-dihydroxy-4-(3-methyl-2-buten-1-yl)phenyl]vinyl}-2,6-dihydroxyphenyl '
                                             'hydrogen sulfate',
                                     'reason': 'No halogen atoms found'}],
    'sample_false_negatives': [   {   'smiles': 'Brn1ccc2ccccc12',
                                      'name': '1-bromoindole',
                                      'reason': 'No carbon-halogen bonds '
                                                'found'},
                                  {   'smiles': 'O=c1cc(oc2c(cccc12)-c1ccccc1)N1CCOCC1',
                                      'name': 'LY294002',
                                      'reason': 'No halogen atoms found'},
                                  {   'smiles': 'Cln1c(=O)n(Cl)c(=O)n(Cl)c1=O',
                                      'name': '1,3,5-trichloro-1,3,5-triazinane-2,4,6-trione',
                                      'reason': 'No carbon-halogen bonds '
                                                'found'},
                                  {   'smiles': '[Zr](Cl)(C1C=CC=C1)(C2C=CC=C2)Cl',
                                      'name': 'Zirconocene dichloride',
                                      'reason': 'No carbon-halogen bonds '
                                                'found'},
                                  {   'smiles': 'O=C(N[C@H](C(=O)O)CC1=CC=CC=C1)[C@H]2N(C(=O)[C@@H](NC(=O)[C@@H](N)[C@H](CC)C)CC=3NC=NC3)CCC2',
                                      'name': 'Angiotensin (5-8)',
                                      'reason': 'No halogen atoms found'},
                                  {   'smiles': 'BrN1C(=O)CCC1=O',
                                      'name': 'N-bromosuccinimide',
                                      'reason': 'No carbon-halogen bonds '
                                                'found'},
                                  {   'smiles': 'CC1=CC(C)=C(NC(=O)C2=CC(=CC(=C2)C#N)S(F)(F)(F)(F)F)C=C1N1C=CN2N=C(C=C12)C1=CC=CN=C1',
                                      'name': 'BAY-826',
                                      'reason': 'No carbon-halogen bonds '
                                                'found'}],
    'attempt': 0,
    'reasoning': 'Let me create a program to classify organohalogen compounds. '
                 'The key requirements are:\n'
                 '\n'
                 '1. Must have at least one carbon-halogen bond\n'
                 '2. Halogens are F, Cl, Br, or I\n'
                 '3. The halogen must be directly bonded to carbon\n'
                 '\n'
                 "Here's the implementation:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 3545,
    'num_false_positives': 100,
    'num_true_negatives': 896,
    'num_false_negatives': 7,
    'num_negatives': None,
    'precision': 0.9725651577503429,
    'recall': 0.9980292792792793,
    'f1': 0.9851326941781298,
    'accuracy': 0.9764731750219877,
    'negative_predictive_value': 0.9922480620155039}