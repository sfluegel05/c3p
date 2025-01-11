"""
Classifies: CHEBI:76983 N-sulfonylurea
"""
"""
Classifies: N-sulfonylurea compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_sulfonylurea(smiles: str):
    """
    Determines if a molecule is an N-sulfonylurea based on its SMILES string.
    N-sulfonylureas contain a urea group where one nitrogen is connected to a sulfonyl group.
    General structure: R-SO2-NH-C(=O)-NH-R'

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-sulfonylurea, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define core N-sulfonylurea pattern
    # [#7] represents any nitrogen
    # The pattern captures R-SO2-NH-C(=O)-NH-R' including variations
    nsulf_pattern = Chem.MolFromSmarts('[#7][#6](=[#8])[#7]S(=O)(=O)[#6,#7,#8]')
    
    # Alternative pattern for reversed arrangement
    nsulf_pattern2 = Chem.MolFromSmarts('[#6,#7,#8]S(=O)(=O)[#7][#6](=[#8])[#7]')
    
    # Pattern for cyclic variations
    cyclic_pattern = Chem.MolFromSmarts('[#7]1[#6](=[#8])[#7]S(=O)(=O)[#6,#7,#8]1')

    if not (mol.HasSubstructMatch(nsulf_pattern) or 
            mol.HasSubstructMatch(nsulf_pattern2) or
            mol.HasSubstructMatch(cyclic_pattern)):
        return False, "Missing required N-sulfonylurea substructure"

    # Exclude invalid cases
    
    # Pattern for direct S-N-N connection without carbonyl
    invalid_pattern1 = Chem.MolFromSmarts('[N]S(=O)(=O)[N]')
    # Pattern for sulfamate derivatives
    invalid_pattern2 = Chem.MolFromSmarts('[O]S(=O)(=O)[N]')
    
    if mol.HasSubstructMatch(invalid_pattern1):
        # Verify that if S-N-N exists, it's part of proper N-sulfonylurea
        for match in mol.GetSubstructMatches(invalid_pattern1):
            n1, s, n2 = match[0], match[1], match[2]
            n1_atom = mol.GetAtomWithIdx(n1)
            n2_atom = mol.GetAtomWithIdx(n2)
            
            # Check if either nitrogen is part of a carbonyl group
            has_carbonyl = False
            for n_atom in [n1_atom, n2_atom]:
                for neighbor in n_atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'C':
                        for c_neighbor in neighbor.GetNeighbors():
                            if c_neighbor.GetSymbol() == 'O' and c_neighbor.GetTotalNumHs() == 0:
                                has_carbonyl = True
                                break
            if not has_carbonyl:
                return False, "Invalid S-N-N arrangement without proper carbonyl group"

    # Additional validation of connectivity
    for match in mol.GetSubstructMatches(nsulf_pattern):
        n1, c, o, n2, s = match[:5]
        # Verify oxidation states and connectivity
        if (mol.GetAtomWithIdx(o).GetTotalNumHs() == 0 and  # Carbonyl oxygen
            mol.GetAtomWithIdx(s).GetTotalValence() == 6):   # Sulfur(VI)
            return True, "Contains N-sulfonylurea group with validated connectivity"
            
    for match in mol.GetSubstructMatches(nsulf_pattern2):
        if len(match) >= 5:  # Ensure we have enough atoms in match
            return True, "Contains N-sulfonylurea group with validated connectivity"
            
    for match in mol.GetSubstructMatches(cyclic_pattern):
        if len(match) >= 5:  # Ensure we have enough atoms in match
            return True, "Contains N-sulfonylurea group with validated connectivity"

    return False, "Structure lacks proper N-sulfonylurea connectivity"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:76983',
                          'name': 'N-sulfonylurea',
                          'definition': 'A urea in which one of the hydrogens '
                                        'attached to a nitrogen of the urea '
                                        'group is replaced by a sulfonyl '
                                        'group. The N-sulfonylurea moiety is a '
                                        'key group in various herbicides, as '
                                        'well as in a number of antidiabetic '
                                        'drugs used in the management of type '
                                        '2 diabetis mellitus.',
                          'parents': ['CHEBI:35358', 'CHEBI:47857'],
                          'xrefs': [   'PMID:23237384',
                                       'PMID:23594109',
                                       'PMID:23609875',
                                       'Wikipedia:Sulfonylurea'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.5528455284552846 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'O=C(NS(CC=1C=CC=CC1C(O)=O)(=O)=O)NC2=NC(=CC(=N2)OC)OC NAME: '
               'bensulfuron REASON: CORRECT Contains N-sulfonylurea group with '
               'validated connectivity\n'
               ' * SMILES: '
               'C1(=NC(=NC(=N1)NC(NS(C=2C(=CC=CC2)C(=O)OC)(=O)=O)=O)OCC)NC '
               'NAME: ethametsulfuron-methyl REASON: CORRECT Contains '
               'N-sulfonylurea group with validated connectivity\n'
               ' * SMILES: Cc1ccc(cc1)S(=O)(=O)NC(=O)NN1CC2CCCC2C1 NAME: '
               'gliclazide REASON: CORRECT Contains N-sulfonylurea group with '
               'validated connectivity\n'
               ' * SMILES: Cc1ccc(cc1)S(=O)(=O)NC(=O)NN1CCCCCC1 NAME: '
               'tolazamide REASON: CORRECT Contains N-sulfonylurea group with '
               'validated connectivity\n'
               ' * SMILES: ClC1=CC=C(S(=O)(=O)NC(=O)NCC(=O)C)C=C1 NAME: '
               '1-(4-chlorobenzenesulfonyl)-3-(2-oxopropyl)urea REASON: '
               'CORRECT Contains N-sulfonylurea group with validated '
               'connectivity\n'
               ' * SMILES: CCCCNC(=O)NS(=O)(=O)c1ccc(C)cc1 NAME: tolbutamide '
               'REASON: CORRECT Contains N-sulfonylurea group with validated '
               'connectivity\n'
               ' * SMILES: '
               'COC1=CC(OC)=NC(NC(=O)NS(=O)(=O)C2=C(C(Cl)=NN2C)C2=NOC[C@H](C)O2)=N1 '
               'NAME: (S)-metazosulfuron REASON: CORRECT Contains '
               'N-sulfonylurea group with validated connectivity\n'
               ' * SMILES: '
               'COC1=CC(OC)=NC(NC(=O)NS(=O)(=O)C2=C(C(Cl)=NN2C)C2=NOCC(C)O2)=N1 '
               'NAME: '
               '3-chloro-N-[(4,6-dimethoxypyrimidin-2-yl)carbamoyl]-1-methyl-4-(5-methyl-5,6-dihydro-1,4,2-dioxazin-3-yl)-1H-pyrazole-5-sulfonamide '
               'REASON: CORRECT Contains N-sulfonylurea group with validated '
               'connectivity\n'
               ' * SMILES: '
               'C1=C(N=C(N=C1OC)NC(NS(C2=C(C=CC=C2)C(OCC)=O)(=O)=O)=O)Cl NAME: '
               'chlorimuron-ethyl REASON: CORRECT Contains N-sulfonylurea '
               'group with validated connectivity\n'
               ' * SMILES: S1C=CC(=C1C(O)=O)S(NC(=O)NC=2N=C(N=C(N2)OC)C)(=O)=O '
               'NAME: thifensulfuron REASON: CORRECT Contains N-sulfonylurea '
               'group with validated connectivity\n'
               ' * SMILES: '
               'O=C(NS(CC=1C=CC=CC1C(OC)=O)(=O)=O)NC2=NC(=CC(=N2)OC)OC NAME: '
               'bensulfuron-methyl REASON: CORRECT Contains N-sulfonylurea '
               'group with validated connectivity\n'
               ' * SMILES: CC1=CC(C)=NC(NC(=O)NS(=O)(=O)C2=CC=CC=C2C(O)=O)=N1 '
               'NAME: sulfometuron REASON: CORRECT Contains N-sulfonylurea '
               'group with validated connectivity\n'
               ' * SMILES: C1(=NC(=NC(=N1)NC(NS(C=2C(=CC=CC2)Cl)(=O)=O)=O)OC)C '
               'NAME: chlorsulfuron REASON: CORRECT Contains N-sulfonylurea '
               'group with validated connectivity\n'
               ' * SMILES: '
               'C1=C(N=C(N=C1OC)NC(NS(C2=C(C=CC=C2)C(O)=O)(=O)=O)=O)Cl NAME: '
               'chlorimuron REASON: CORRECT Contains N-sulfonylurea group with '
               'validated connectivity\n'
               ' * SMILES: '
               'C1=CC=C(C(=C1C(OC)=O)S(NC(NC2=NC(=NC(=N2)N(C)C)OCC(F)(F)F)=O)(=O)=O)C '
               'NAME: triflusulfuron-methyl REASON: CORRECT Contains '
               'N-sulfonylurea group with validated connectivity\n'
               ' * SMILES: '
               'COC1=CC(OC)=NC(NC(=O)NS(=O)(=O)C2=C(C(Cl)=NN2C)C2=NOC[C@@H](C)O2)=N1 '
               'NAME: (R)-metazosulfuron REASON: CORRECT Contains '
               'N-sulfonylurea group with validated connectivity\n'
               ' * SMILES: CCCCCN(C(=O)NC(=O)NC1C2CC3CC1CC(C2)C3)S(=O)(=O)C '
               'NAME: '
               '3-[(2-adamantylamino)-oxomethyl]-1-methylsulfonyl-1-pentylurea '
               'REASON: CORRECT Contains N-sulfonylurea group with validated '
               'connectivity\n'
               ' * SMILES: COC(=O)c1ccccc1S(=O)(=O)NC(=O)Nc1nc(C)nc(OC)n1 '
               'NAME: metsulfuron methyl REASON: CORRECT Contains '
               'N-sulfonylurea group with validated connectivity\n'
               ' * SMILES: '
               'CC(C)(O)C1=COC(=C1)S(=O)(=O)NC(=O)NC1=C2CCCC2=CC2=C1CCC2 NAME: '
               'MCC950 REASON: CORRECT Contains N-sulfonylurea group with '
               'validated connectivity\n'
               ' * SMILES: '
               'C1(=NC(C)=NC(OC)=N1)N(C(=O)NS(C=2C(C(O)=O)=CC=CC2)(=O)=O)C '
               'NAME: tribenuron REASON: CORRECT Contains N-sulfonylurea group '
               'with validated connectivity\n'
               ' * SMILES: '
               'C=1(C(NC=2C=C(C=CC2)C(O)=O)=CC=NC1)S(NC(NC(C)C)=O)(=O)=O NAME: '
               'torasemide carboxylic acid REASON: CORRECT Contains '
               'N-sulfonylurea group with validated connectivity\n'
               ' * SMILES: COCCOc1ccccc1S(=O)(=O)NC(=O)Nc1nc(OC)nc(OC)n1 NAME: '
               'Cinosulfuron REASON: CORRECT Contains N-sulfonylurea group '
               'with validated connectivity\n'
               ' * SMILES: '
               'Cc1cnc(cn1)C(=O)NCCc1ccc(cc1)S(=O)(=O)NC(=O)NC1CCCCC1 NAME: '
               'glipizide REASON: CORRECT Contains N-sulfonylurea group with '
               'validated connectivity\n'
               ' * SMILES: '
               'COc1ccc(Cl)cc1C(=O)NCCc1ccc(cc1)S(=O)(=O)NC(=O)NC1CCCCC1 NAME: '
               'glyburide REASON: CORRECT Contains N-sulfonylurea group with '
               'validated connectivity\n'
               ' * SMILES: '
               'S1C=CC(=C1C(OC)=O)S(NC(=O)NC=2N=C(N=C(N2)OC)C)(=O)=O NAME: '
               'thifensulfuron-methyl REASON: CORRECT Contains N-sulfonylurea '
               'group with validated connectivity\n'
               'False positives: SMILES: '
               'COc1cc(OC)nc(NC(=O)NS(=O)(=O)N(C)S(C)(=O)=O)n1 NAME: '
               'Amidosulfuron REASON: WRONGLY CLASSIFIED Contains '
               'N-sulfonylurea group with validated connectivity\n'
               ' * SMILES: S(NC(NC1CCCCC1)=O)(=O)(=O)C2=CC3=C(C=C2)CCC3 NAME: '
               'glyhexamide REASON: WRONGLY CLASSIFIED Contains N-sulfonylurea '
               'group with validated connectivity\n'
               ' * SMILES: S(=O)(=O)(NC(=O)NN1CC2C(CC(O)C2)C1)C3=CC=C(C=C3)C '
               'NAME: '
               'N-((5-hydroxyhexahydrocyclopenta[c]pyrrol-2(1H)-yl)carbamoyl)-4-methylbenzenesulfonamide '
               'REASON: WRONGLY CLASSIFIED Contains N-sulfonylurea group with '
               'validated connectivity\n'
               ' * SMILES: '
               'COc1cc(OC)nc(NC(=O)NS(=O)(=O)c2cc(NC=O)ccc2C(=O)N(C)C)n1 NAME: '
               'foramsulfuron REASON: WRONGLY CLASSIFIED Contains '
               'N-sulfonylurea group with validated connectivity\n'
               ' * SMILES: '
               'S(=O)(=O)(NC(=O)N[C@@H]1C[C@H](O)CCC1)C2=CC=C(CCNC(=O)C3=NC=C(N=C3)C)C=C2 '
               'NAME: 3-cis-Hydroxyglipizide REASON: WRONGLY CLASSIFIED '
               'Contains N-sulfonylurea group with validated connectivity\n'
               ' * SMILES: C1=C(C=CC(=C1)S(NC(NN2CCCC2)=O)(=O)=O)Cl NAME: '
               'Glyclopyramide REASON: WRONGLY CLASSIFIED Contains '
               'N-sulfonylurea group with validated connectivity\n'
               ' * SMILES: S(=O)(=O)(NC(=O)N1CCCC1)C2=CC=C(C=C2)C NAME: '
               'tolpyrramide REASON: WRONGLY CLASSIFIED Contains '
               'N-sulfonylurea group with validated connectivity\n'
               ' * SMILES: '
               'S(=O)(=O)(NC(=O)NCCCCC)C1=C(NC2CCCCC2)C=CC([N+]([O-])=O)=C1 '
               'NAME: BM 567 REASON: WRONGLY CLASSIFIED Contains '
               'N-sulfonylurea group with validated connectivity\n'
               ' * SMILES: S(NC(NC1CCCCC1)=O)(=O)(=O)C2=CC(N)=C(C=C2)C NAME: '
               'metahexamide REASON: WRONGLY CLASSIFIED Contains '
               'N-sulfonylurea group with validated connectivity\n'
               ' * SMILES: S(=O)(=O)(NC(=O)NC1CCCCC1)C2=CC=C(C(O)C)C=C2 NAME: '
               'Hydroxyhexamide REASON: WRONGLY CLASSIFIED Contains '
               'N-sulfonylurea group with validated connectivity\n'
               ' * SMILES: ClC1=CC=C(S(=O)(=O)NC(=O)NCCCO)C=C1 NAME: '
               '1-(4-chlorophenyl)sulfonyl-3-(3-hydroxypropyl)urea REASON: '
               'WRONGLY CLASSIFIED Contains N-sulfonylurea group with '
               'validated connectivity\n'
               ' * SMILES: CCC1=CC(=NC=C1)NC(=O)NS(=O)(=O)C2=CC=C(C=C2)C NAME: '
               '1-(4-ethyl-2-pyridinyl)-3-(4-methylphenyl)sulfonylurea REASON: '
               'WRONGLY CLASSIFIED Contains N-sulfonylurea group with '
               'validated connectivity\n'
               ' * SMILES: COC1=NC(C)=NC(NC(=O)[N-]S(=O)(=O)C2=CC=CC=C2I)=N1 '
               'NAME: iofensulfuron(1-) REASON: WRONGLY CLASSIFIED Contains '
               'N-sulfonylurea group with validated connectivity\n'
               ' * SMILES: S(=O)(=O)(NC(=O)N1CCOCC1)C2=CC=C(C=C2)C NAME: '
               'N1-morpholinocarbonyl-4-methylbenzene-1-sulfonamide REASON: '
               'WRONGLY CLASSIFIED Contains N-sulfonylurea group with '
               'validated connectivity\n'
               ' * SMILES: '
               '[H][C@]12SC(C)(C)[C@@H](N1C(=O)[C@H]2NC(=O)[C@H](NC(=O)N1CCN(C1=O)S(C)(=O)=O)c1ccccc1)C([O-])=O '
               'NAME: mezlocillin(1-) REASON: WRONGLY CLASSIFIED Contains '
               'N-sulfonylurea group with validated connectivity\n'
               ' * SMILES: '
               'S(=O)(=O)(NC(=O)NC1CCC(CC1)CO)C2=CC=C(CCNC(=O)N3CC(=C(C3=O)CC)C)C=C2 '
               'NAME: trans-Hydroxy Glimepiride REASON: WRONGLY CLASSIFIED '
               'Contains N-sulfonylurea group with validated connectivity\n'
               ' * SMILES: COc1cc(OC)nc(NC(=O)NS(=O)(=O)c2c(Cl)nc3ccccn23)n1 '
               'NAME: Imazosulfuron REASON: WRONGLY CLASSIFIED Contains '
               'N-sulfonylurea group with validated connectivity\n'
               ' * SMILES: '
               'ClC1=CC(=C(OC)C=C1)C(=O)NCCC2=CC=C(S(=O)(=O)NC(=O)NC3CCC(O)CC3)C=C2 '
               'NAME: trans-4-Hydroxyglyburide REASON: WRONGLY CLASSIFIED '
               'Contains N-sulfonylurea group with validated connectivity\n'
               ' * SMILES: CCOc1ccccc1OS(=O)(=O)NC(=O)Nc1nc(OC)cc(OC)n1 NAME: '
               'Ethoxysulfuron REASON: WRONGLY CLASSIFIED Contains '
               'N-sulfonylurea group with validated connectivity\n'
               ' * SMILES: CCCCNC(=O)NS(=O)(=O)c1ccc(CO)cc1 NAME: '
               '4-hydroxytolbutamide REASON: WRONGLY CLASSIFIED Contains '
               'N-sulfonylurea group with validated connectivity\n'
               ' * SMILES: '
               'C[C@@]12C(C)(C)[C@@H]([C@@H]([C@@H]1O)NC(NS(=O)(=O)C3=CC=C(C=C3)C)=O)CC2 '
               'NAME: glibornuride REASON: WRONGLY CLASSIFIED Contains '
               'N-sulfonylurea group with validated connectivity\n'
               ' * SMILES: '
               '[NaH].COC(=O)c1ccc(nc1S(=O)(=O)NC(=O)Nc1nc(OC)cc(OC)n1)C(F)(F)F '
               'NAME: Flupyrsulfuron-methyl sodium REASON: WRONGLY CLASSIFIED '
               'Contains N-sulfonylurea group with validated connectivity\n'
               ' * SMILES: CCCCNC(=O)NS(=O)(=O)C1=CC=C(N)C=C1 NAME: '
               'carbutamide REASON: WRONGLY CLASSIFIED Contains N-sulfonylurea '
               'group with validated connectivity\n'
               ' * SMILES: '
               '[Na+].[H][C@]12SC(C)(C)[C@@H](N1C(=O)[C@H]2NC(=O)[C@H](NC(=O)N1CCN(C1=O)S(C)(=O)=O)c1ccccc1)C([O-])=O '
               'NAME: mezlocillin sodium REASON: WRONGLY CLASSIFIED Contains '
               'N-sulfonylurea group with validated connectivity\n'
               ' * SMILES: S(=O)(=O)(NC(=O)N1CCN(CC1)C)C2=CC=C(C=C2)C NAME: '
               '4-Methyl-N-((4-methyl-1-piperazinyl)carbonyl)benzenesulfonamide '
               'REASON: WRONGLY CLASSIFIED Contains N-sulfonylurea group with '
               'validated connectivity\n'
               'False negatives: SMILES: '
               'C=1(N(C(N(N1)C(NS(C2=CC=CC=C2OC(F)(F)F)(=O)=O)=O)=O)C)OC NAME: '
               'flucarbazone REASON: MISSED Missing required N-sulfonylurea '
               'substructure\n'
               ' * SMILES: '
               'S(NC(=O)N1C(N(C(=N1)OC)C)=O)(C2=C(SC=C2C(=O)O)C)(=O)=O NAME: '
               'thiencarbazone REASON: MISSED Missing required N-sulfonylurea '
               'substructure\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'CCS(=O)(=O)NC1=CC2=C(C=C1)OC[C@H](N(C[C@H]([C@H](CN(C2=O)C)OC)C)C(=O)COC)C',
                                     'name': 'N-[(4R,7R,8R)-8-methoxy-5-(2-methoxy-1-oxoethyl)-4,7,10-trimethyl-11-oxo-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]ethanesulfonamide',
                                     'reason': 'Missing required '
                                               'N-sulfonylurea substructure'},
                                 {   'smiles': 'C1CCC(C1)(C2=CN(N=N2)CC[C@@H]3CC[C@@H]([C@H](O3)CO)NC(=O)C4=CC=CC=C4F)O',
                                     'name': '2-fluoro-N-[(2S,3S,6S)-6-[2-[4-(1-hydroxycyclopentyl)-1-triazolyl]ethyl]-2-(hydroxymethyl)-3-oxanyl]benzamide',
                                     'reason': 'Missing required '
                                               'N-sulfonylurea substructure'},
                                 {   'smiles': 'S(=O)(=O)(NC)CC1=CC2=C(NC=C2CC(O)=O)C=C1',
                                     'name': '1H-Indole-3-acetic acid, '
                                             '5-[[(methylamino)sulfonyl]methyl]-',
                                     'reason': 'Missing required '
                                               'N-sulfonylurea substructure'},
                                 {   'smiles': 'O([C@@H]1O[C@@H]([C@H](O)[C@H](O[C@H]2O[C@@H]([C@H](O)[C@H](O)[C@H]2NC(=O)C)CO)[C@H]1O[C@@H]3O[C@H]([C@@H](O)[C@@H](O)[C@@H]3O)C)CO)[C@H]4[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]4CO)OC[C@H]5O[C@H](O)[C@H](NC(=O)C)[C@@H](O[C@@H]6O[C@@H]([C@H](O)[C@H](O)[C@H]6O[C@@H]7O[C@H]([C@@H](O)[C@@H](O)[C@@H]7O)C)CO)[C@H]5O',
                                     'name': 'N-[(2R,3R,4R,5R,6R)-2-[(2S,3R,4S,5S,6R)-2-[(2R,3S,4R,5R,6R)-5-Acetamido-6-[[(2R,3R,4R,5R,6S)-5-acetamido-4-[(2R,3R,4S,5R,6R)-4,5-dihydroxy-6-(hydroxymethyl)-3-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-2-yl]oxy-3,6-dihydroxyoxan-2-yl]methoxy]-4-hydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-5-hydroxy-6-(hydroxymethyl)-3-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-4-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]acetamide',
                                     'reason': 'Missing required '
                                               'N-sulfonylurea substructure'},
                                 {   'smiles': 'O([C@H]1[C@H](O)[C@H](O[C@@H](O[C@H]2[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]2CO)O[C@@H]([C@H](O)[C@@H](NC(=O)C)CO)[C@H](O)CO)[C@H]1O)CO[C@H]3O[C@@H]([C@@H](O)[C@H](O[C@H]4O[C@@H]([C@@H](O)[C@H](O)[C@@H]4O[C@H]5O[C@@H]([C@@H](O)[C@H](O)[C@@H]5O)CO)CO)[C@@H]3O)CO[C@H]6O[C@@H]([C@@H](O)[C@H](O)[C@@H]6O[C@H]7O[C@@H]([C@@H](O)[C@H](O)[C@@H]7O)CO)CO)[C@H]8O[C@@H]([C@@H](O)[C@H](O)[C@@H]8O[C@H]9O[C@@H]([C@@H](O)[C@H](O)[C@@H]9O[C@H]%10O[C@@H]([C@@H](O)[C@H](O)[C@@H]%10O)CO)CO)CO',
                                     'name': 'N-[(2S,3R,4R,5S,6R)-2-[(2R,3S,4R,5S)-5-Acetamido-1,2,4,6-tetrahydroxyhexan-3-yl]oxy-5-[(2S,3S,4S,5R,6R)-4-[(2R,3S,4S,5S,6R)-3-[(2S,3S,4S,5S,6R)-4,5-dihydroxy-6-(hydroxymethyl)-3-[(2R,3S,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-6-[[(2S,3S,4S,5R,6R)-4-[(2S,3S,4S,5S,6R)-4,5-dihydroxy-6-(hydroxymethyl)-3-[(2R,3S,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxy-6-[[(2S,3S,4S,5S,6R)-4,5-dihydroxy-6-(hydroxymethyl)-3-[(2R,3S,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxymethyl]-3,5-dihydroxyoxan-2-yl]oxymethyl]-3,5-dihydroxyoxan-2-yl]oxy-4-hydroxy-6-(hydroxymethyl)oxan-3-yl]acetamide',
                                     'reason': 'Missing required '
                                               'N-sulfonylurea substructure'},
                                 {   'smiles': 'P(=O)(O[C@H]([C@](O)(C=C[C@@H]1OC(=O)C=C[C@@H]1CC)CCN)C[C@H](O)C=CC=C[C@H]2C[C@@H](OC(=O)C3CCCCC3)CCC2)(O)O',
                                     'name': 'Phoslactomycin E',
                                     'reason': 'Missing required '
                                               'N-sulfonylurea substructure'},
                                 {   'smiles': 'ClC1=C(Cl)C=CC=C1C2=N[N+](O)=C(N)N=C2N',
                                     'name': '3,5-diamino-6-(2,3-dichlorophenyl)-2-hydroxy-1,2,4-triazin-2-ium',
                                     'reason': 'Missing required '
                                               'N-sulfonylurea substructure'},
                                 {   'smiles': 'P(CC[N+](C)(C)C)(OCC(OC(=O)C)COC(=O)C)(=O)[O-]',
                                     'name': 'diacylglycerol '
                                             '2-trimethylaminoethylphosphonate',
                                     'reason': 'Missing required '
                                               'N-sulfonylurea substructure'},
                                 {   'smiles': 'C[C@@H]1CN([C@@H](COC2=C(C=CC(=C2)NC(=O)COC)C(=O)N(C[C@@H]1OC)C)C)CC3=CC(=CC=C3)F',
                                     'name': 'N-[(5R,6R,9R)-8-[(3-fluorophenyl)methyl]-5-methoxy-3,6,9-trimethyl-2-oxo-11-oxa-3,8-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]-2-methoxyacetamide',
                                     'reason': 'Missing required '
                                               'N-sulfonylurea substructure'},
                                 {   'smiles': 'NCC(=O)N1[C@H](C(=O)O)C[C@H](C1)O',
                                     'name': 'Gly-Hyp',
                                     'reason': 'Missing required '
                                               'N-sulfonylurea substructure'}],
    'sample_false_negatives': [],
    'attempt': 2,
    'reasoning': 'Let me analyze the previous outcomes and suggest '
                 'improvements:\n'
                 '\n'
                 '1. Pattern Issues:\n'
                 '- The program had too many false positives, suggesting the '
                 'patterns were too permissive\n'
                 '- The N-sulfonylurea structure needs to be more specifically '
                 'defined as -N(H)-C(=O)-N(H)-S(=O)(=O)-R\n'
                 "- The current patterns allow for variations that shouldn't "
                 'be included\n'
                 '\n'
                 '2. Connectivity Validation:\n'
                 '- The current validation checks for N-S and N-C connectivity '
                 "but doesn't properly ensure the specific arrangement\n"
                 '- Need to verify the exact sequence of atoms in the '
                 'N-sulfonylurea group\n'
                 '- Should check that one nitrogen is connected to both '
                 'carbonyl and sulfonyl groups\n'
                 '\n'
                 '3. False Negatives Analysis:\n'
                 '- Some valid compounds (flucarbazone, thiencarbazone) were '
                 'missed\n'
                 '- These compounds have slightly different arrangements but '
                 'are still N-sulfonylureas\n'
                 '- Need to account for cyclic variations of the basic '
                 'structure\n'
                 '\n'
                 "Here's my improved version:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 36,
    'num_false_positives': 54,
    'num_true_negatives': 142210,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.4,
    'recall': 1.0,
    'f1': 0.5714285714285715,
    'accuracy': 0.9996205200281096,
    'negative_predictive_value': 1.0}