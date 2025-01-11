"""
Classifies: CHEBI:16385 organic sulfide
"""
"""
Classifies: organic sulfide compounds (thioethers)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_organic_sulfide(smiles: str):
    """
    Determines if a molecule is an organic sulfide based on its SMILES string.
    An organic sulfide has the structure R-S-R where R is any organic group (not H).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organic sulfide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find sulfur atoms
    sulfur_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16]
    if not sulfur_atoms:
        return False, "No sulfur atoms found"

    # Check each sulfur atom
    valid_sulfide_count = 0
    for sulfur in sulfur_atoms:
        # Get neighboring atoms
        neighbors = sulfur.GetNeighbors()
        
        # Count number of carbon neighbors
        carbon_neighbors = sum(1 for n in neighbors if n.GetAtomicNum() == 6)
        
        # Check oxidation state (avoid sulfoxides and sulfones)
        explicit_valence = sulfur.GetExplicitValence()
        
        # Check if sulfur has exactly two carbon neighbors and proper valence
        if carbon_neighbors == 2 and explicit_valence == 2:
            valid_sulfide_count += 1
            
        # Check for special cases like thiocyanates (-SCN)
        elif len(neighbors) == 2 and any(n.GetAtomicNum() == 6 for n in neighbors):
            if any(n.GetAtomicNum() == 7 and any(nb.GetAtomicNum() == 6 for nb in n.GetNeighbors()) 
                  for n in neighbors):
                valid_sulfide_count += 1

    if valid_sulfide_count == 0:
        # Check specific patterns that might have been missed
        sulfide_pattern = Chem.MolFromSmarts("[#6]-[#16X2]-[#6]")
        if mol.HasSubstructMatch(sulfide_pattern):
            return True, "Contains R-S-R pattern where R is carbon-containing group"
        return False, "No sulfur atom with two carbon neighbors in correct oxidation state"

    return True, f"Contains {valid_sulfide_count} sulfide group(s) with R-S-R pattern"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:16385',
                          'name': 'organic sulfide',
                          'definition': 'Compounds having the structure RSR (R '
                                        '=/= H). Such compounds were once '
                                        'called thioethers.',
                          'parents': ['CHEBI:26822', 'CHEBI:33261'],
                          'xrefs': ['KEGG:C00297'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'S(=O)(=O)(C1=C2C(=CC3=C1NC=4C=CC=CC34)[C@@]5([C@H]([C@](C(=O)O)([C@@H](O)CC5)C)CC2)C)C6=CC7=C(NC8=C7C=C9[C@@]%10([C@H]([C@](C(=O)O)([C@@H](O)CC%10)C)CCC9=C8)C)C=C6',
                                     'name': 'Sulfadixiamycin C',
                                     'reason': 'No sulfur atom with two carbon '
                                               'neighbors in correct oxidation '
                                               'state'},
                                 {   'smiles': 'CNC(O)=O',
                                     'name': 'methylcarbamic acid',
                                     'reason': 'No sulfur atoms found'},
                                 {   'smiles': 'CCNC(=O)NC1=CC2=C(C=C1)OC[C@H]3[C@@H](CC[C@H](O3)CC(=O)N[C@@H](C)C4=CC=CC=C4)N(C2=O)C',
                                     'name': '2-[(2S,4aR,12aR)-8-(ethylcarbamoylamino)-5-methyl-6-oxo-2,3,4,4a,12,12a-hexahydropyrano[2,3-c][1,5]benzoxazocin-2-yl]-N-[(1S)-1-phenylethyl]acetamide',
                                     'reason': 'No sulfur atoms found'},
                                 {   'smiles': 'O([C@H]1O[C@@H]([C@@H](O[C@@H]2O[C@@H]([C@@H](O[C@@H]3O[C@@H]([C@H](O)[C@H](O)[C@H]3O)CO[C@]4(O[C@H]([C@H](NC(=O)C)[C@@H](O)C4)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]2NC(=O)C)CO)[C@H](O)[C@@H]1O[C@@H]5O[C@@H]([C@@H](O[C@@H]6O[C@@H]([C@H](O)[C@H](O)[C@H]6O)CO[C@]7(O[C@H]([C@H](NC(=O)C)[C@@H](O)C7)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]5NC(=O)C)CO)CO)[C@H]8[C@H](O)[C@H](O[C@@H](O[C@H]9[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]9CO)O[C@H]%10[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]%10CO[C@@H]%11O[C@H]([C@@H](O)[C@@H](O)[C@@H]%11O)C)O)[C@H]8O)CO[C@H]%12O[C@@H]([C@@H](O)[C@H](O)[C@@H]%12O[C@@H]%13O[C@@H]([C@@H](O[C@@H]%14O[C@@H]([C@H](O)[C@H](O[C@]%15(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%15)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]%14O)CO)[C@H](O)[C@H]%13NC(=O)C)CO)CO[C@@H]%16O[C@@H]([C@@H](O[C@@H]%17O[C@@H]([C@H](O)[C@H](O)[C@H]%17O)CO[C@]%18(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%18)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]%16NC(=O)C)CO',
                                     'name': 'CID 91851985',
                                     'reason': 'No sulfur atoms found'},
                                 {   'smiles': 'O(C(=O)C(C1C(CN2C(C1)C=3NC=4C(C3CC2)=CC=CC4)CC)=COC)C',
                                     'name': 'Methyl '
                                             '2-(3-ethyl-1,2,3,4,6,7,12,12b-octahydroindolo[2,3-a]quinolizin-2-yl)-3-methoxyprop-2-enoate',
                                     'reason': 'No sulfur atoms found'},
                                 {   'smiles': 'O[C@H](/C=C/C=C/C=C/[C@H](O)[C@H](O)C=C)[C@H](O)/C=C/C',
                                     'name': 'Separacene C',
                                     'reason': 'No sulfur atoms found'},
                                 {   'smiles': 'C[C@@H]1O[C@@H](O[C@@H]2[C@@H](CO)O[C@@H](O[C@@H]3[C@@H](O)C(O)O[C@H](CO)[C@@H]3O)[C@H](NC(C)=O)[C@H]2O[C@@H]2O[C@H](CO)[C@H](O)[C@H](OS(O)(=O)=O)[C@H]2O)[C@@H](O)[C@H](O)[C@@H]1O',
                                     'name': 'alpha-L-Fucp-(1->4)-[beta-D-Galp3S-(1->3)]-beta-D-GlcpNAc-(1->3)-D-Galp',
                                     'reason': 'No sulfur atom with two carbon '
                                               'neighbors in correct oxidation '
                                               'state'},
                                 {   'smiles': 'C1=CC=CC2=C1C(N([C@H](C(N2)=O)CC=3C=CC(=CC3)OC)C)=O',
                                     'name': "(S)-4'-methoxycyclopeptine",
                                     'reason': 'No sulfur atoms found'},
                                 {   'smiles': 'O=C(N[C@@H](CC=1C=2C(NC1)=CC=CC2)C(O)=O)[C@@H](NC(=O)[C@@H](N)C(C)C)C(C)C',
                                     'name': 'Val-Val-Trp',
                                     'reason': 'No sulfur atoms found'},
                                 {   'smiles': 'C=1C(=C(C=CC1/C=C/CO)OC(CO)C(O)C=2C=C(C(=CC2)O)OC)OC',
                                     'name': 'guaiacylglycerol beta-coniferyl '
                                             'ether',
                                     'reason': 'No sulfur atoms found'}],
    'sample_false_negatives': [   {   'smiles': 'C1=2C(=CNC1=CC=CC2)C\\C(\\[S-])=N\\OS([O-])(=O)=O',
                                      'name': '(Z)-(indol-3-yl)-N-(sulfonatooxy)methanimidothioate',
                                      'reason': 'No sulfur atom with two '
                                                'carbon neighbors in correct '
                                                'oxidation state'},
                                  {   'smiles': 'C1=CC=C(C=C1)CCCNC2=[NH+]C=NC3=C2NC=N3',
                                      'name': 'LSM-21130',
                                      'reason': 'No sulfur atoms found'},
                                  {   'smiles': 'CCCS(=O)(=O)c1ccc2[nH]c(NC(=O)OC)nc2c1',
                                      'name': 'Albendazole sulfone',
                                      'reason': 'No sulfur atom with two '
                                                'carbon neighbors in correct '
                                                'oxidation state'},
                                  {   'smiles': 'CCC(=O)N(C)C1=CC=C(C=C1)NC2=NC3=C(C(=N2)NC4CCCN(C4)C(=O)C=C)NC=N3',
                                      'name': 'N-methyl-N-[4-[[6-[[1-(1-oxoprop-2-enyl)-3-piperidinyl]amino]-7H-purin-2-yl]amino]phenyl]propanamide',
                                      'reason': 'No sulfur atoms found'},
                                  {   'smiles': '[S-]/C(=N\\OS([O-])(=O)=O)/CC=C',
                                      'name': '(Z)-N-(sulfonatooxy)prop-2-enimidothioate',
                                      'reason': 'No sulfur atom with two '
                                                'carbon neighbors in correct '
                                                'oxidation state'},
                                  {   'smiles': 'O=C1NC=2N=CN=C(C2N1)NCC=C(C)C',
                                      'name': '8-oxo-isopentenyladenine',
                                      'reason': 'No sulfur atoms found'},
                                  {   'smiles': 'N1(C2=C(C(=NC=N2)S)N=C1)[C@@H]3O[C@H](COP(O)(O)=O)[C@H]([C@H]3O)O',
                                      'name': '9-(5-phosphoribofuranosyl)-6-mercaptopurine',
                                      'reason': 'No sulfur atom with two '
                                                'carbon neighbors in correct '
                                                'oxidation state'},
                                  {   'smiles': 'SC(S)=N',
                                      'name': 'carbonimidodithioic acid',
                                      'reason': 'No sulfur atom with two '
                                                'carbon neighbors in correct '
                                                'oxidation state'},
                                  {   'smiles': 'CC(C)N1C=NC2=C1N=C(N=C2NCC3=CC=CC=C3)NCCO',
                                      'name': '2-[[6-[(phenylmethyl)amino]-9-propan-2-yl-2-purinyl]amino]ethanol',
                                      'reason': 'No sulfur atoms found'},
                                  {   'smiles': 'COC1=CC=CC=C1CNC2=[NH+]C=NC3=C2NC=N3',
                                      'name': 'LSM-32859',
                                      'reason': 'No sulfur atoms found'}],
    'attempt': 0,
    'reasoning': 'Let me help create a program to classify organic sulfides. '
                 'Based on the definition and examples, here are the key '
                 'characteristics:\n'
                 '\n'
                 '1. Must contain at least one sulfur atom bonded to two '
                 'carbon atoms (R-S-R where R ≠ H)\n'
                 '2. Should not count thiols (R-SH), sulfoxides (R-S(=O)-R), '
                 'or sulfones (R-S(=O)2-R)\n'
                 '3. Can include aromatic and aliphatic carbon groups\n'
                 '4. May contain other functional groups and heteroatoms\n'
                 '\n'
                 "Here's the implementation:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 564,
    'num_false_positives': 100,
    'num_true_negatives': 1589,
    'num_false_negatives': 25,
    'num_negatives': None,
    'precision': 0.8493975903614458,
    'recall': 0.9575551782682513,
    'f1': 0.9002394253790903,
    'accuracy': 0.9451273046532046,
    'negative_predictive_value': 0.9845105328376704}