"""
Classifies: CHEBI:13248 anilide
"""
"""
Classifies: anilide compounds
An anilide is an aromatic amide obtained by acylation of aniline.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_anilide(smiles: str):
    """
    Determines if a molecule is an anilide based on its SMILES string.
    An anilide has a phenyl ring connected to an amide nitrogen.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an anilide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Look for N-phenyl amide pattern
    # [c] - aromatic carbon (part of phenyl ring)
    # [NX3] - nitrogen with 3 connections
    # [CX3](=[OX1]) - amide carbon with double bonded oxygen
    anilide_pattern = Chem.MolFromSmarts('[c][NX3][CX3](=[OX1])')
    
    matches = mol.GetSubstructMatches(anilide_pattern)
    if not matches:
        return False, "No N-phenyl amide group found"
        
    # For each match, verify the aromatic ring is a proper phenyl or substituted phenyl
    for match in matches:
        arom_c = match[0]  # First atom in match is the aromatic carbon
        ring = Chem.GetSymmSSSR(mol)
        
        # Check each ring that contains our aromatic carbon
        for r in ring:
            if arom_c in r:
                # Get all atoms in the ring
                ring_atoms = set(r)
                ring_size = len(ring_atoms)
                
                # Count aromatic atoms
                aromatic_count = sum(1 for atom_idx in ring_atoms 
                                   if mol.GetAtomWithIdx(atom_idx).GetIsAromatic())
                
                # Verify it's a 6-membered aromatic ring
                if ring_size == 6 and aromatic_count == 6:
                    # Verify ring atoms are C or substituted C
                    is_phenyl = True
                    for atom_idx in ring_atoms:
                        atom = mol.GetAtomWithIdx(atom_idx)
                        if atom.GetAtomicNum() != 6:  # not carbon
                            is_phenyl = False
                            break
                    
                    if is_phenyl:
                        return True, "Contains N-phenyl amide group"
                        
    return False, "No proper phenyl ring attached to amide nitrogen"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:13248',
                          'name': 'anilide',
                          'definition': 'Any aromatic amide obtained by '
                                        'acylation of aniline.',
                          'parents': ['CHEBI:22712', 'CHEBI:62733'],
                          'xrefs': [   'KEGG:C01402',
                                       'PMID:23535982',
                                       'PMID:23968552',
                                       'PMID:24273122',
                                       'PMID:6205897'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.2857142857142857 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'FC=1C(N(C2CCN(CC2)CCC3=CC=CC=C3)C(=O)COC)=CC=CC1 NAME: '
               'Ocfentanil REASON: CORRECT Contains N-phenyl amide group\n'
               ' * SMILES: CCC(=O)N(c1ccccc1)C1(CCN(CCc2cccs2)CC1)COC NAME: '
               'sufentanil REASON: CORRECT Contains N-phenyl amide group\n'
               ' * SMILES: CCC(=O)N(c1ccccc1)C1(CCN(CCC(=O)OC)CC1)C(=O)OC '
               'NAME: remifentanil REASON: CORRECT Contains N-phenyl amide '
               'group\n'
               ' * SMILES: CCNC1=C(N=C(C(=N1)C)CN(C2=CC=CC=C2)C(=O)C)C#N NAME: '
               'N-[[6-cyano-5-(ethylamino)-3-methyl-2-pyrazinyl]methyl]-N-phenylacetamide '
               'REASON: CORRECT Contains N-phenyl amide group\n'
               ' * SMILES: '
               'O=C(N(C(C1=CC=C(OC)C=C1)C(=O)NCC2OCCC2)C3=CC=C(NC(=O)C)C=C3)CN4N=NC=5C4=CC=CC5 '
               'NAME: CID 3191469 REASON: CORRECT Contains N-phenyl amide '
               'group\n'
               ' * SMILES: CCC(=O)N(C1CCN(CCc2cccs2)CC1C)c1ccccc1 NAME: '
               '3-methylthiofentanyl REASON: CORRECT Contains N-phenyl amide '
               'group\n'
               ' * SMILES: '
               'C(CC(O)=O)(CC(O)=O)(C(O)=O)O.C1=CC=CC(=C1)N(C2(CCN(CC2)CCC3=CC=CS3)COC)C(=O)CC '
               'NAME: Sufentanil citrate REASON: CORRECT Contains N-phenyl '
               'amide group\n'
               ' * SMILES: '
               'CN(C1=CC=CC=C1)C(=O)COC2=NC(=C(N=N2)C3=CC=CS3)C4=CC=CS4 NAME: '
               '2-[(5,6-dithiophen-2-yl-1,2,4-triazin-3-yl)oxy]-N-methyl-N-phenylacetamide '
               'REASON: CORRECT Contains N-phenyl amide group\n'
               ' * SMILES: CCCOCCN(C(=O)CCl)c1c(CC)cccc1CC NAME: Pretilachlor '
               'REASON: CORRECT Contains N-phenyl amide group\n'
               ' * SMILES: CCC(=O)N(C1CCNCC1)c1ccccc1 NAME: norfentanyl '
               'REASON: CORRECT Contains N-phenyl amide group\n'
               ' * SMILES: N(C(CCCC(O)=O)=O)(C)C1=C(I)C(N)=C(C=C1I)I NAME: '
               'iomeglamic acid REASON: CORRECT Contains N-phenyl amide group\n'
               ' * SMILES: CC(C)N(C(=O)CCl)c1ccccc1 NAME: propachlor REASON: '
               'CORRECT Contains N-phenyl amide group\n'
               ' * SMILES: CCC(=O)N(C1CCN(CC1)CCc1ccccc1)c1ccccc1 NAME: '
               'fentanyl REASON: CORRECT Contains N-phenyl amide group\n'
               ' * SMILES: '
               'N(C(COC1=CC=C(C=C1)OCCCC)=O)(CCN(CC)CC)C2=C(OCC)C=CC(=C2)OCC '
               'NAME: fenoxedil REASON: CORRECT Contains N-phenyl amide group\n'
               ' * SMILES: CC1=CC(=CC=C1)NC(=O)CC2C(=O)N(C(=O)S2)C3=CC=CC=C3 '
               'NAME: '
               '2-(2,4-dioxo-3-phenyl-5-thiazolidinyl)-N-(3-methylphenyl)acetamide '
               'REASON: CORRECT Contains N-phenyl amide group\n'
               ' * SMILES: CCCCCCCN(C1=CC=CC(=C1)C)C(=O)CCl NAME: '
               '2-chloro-N-heptyl-N-(3-methylphenyl)acetamide REASON: CORRECT '
               'Contains N-phenyl amide group\n'
               ' * SMILES: CC(C)COCN(C(=O)CCl)c1c(C)cccc1C NAME: Delachlor '
               'REASON: CORRECT Contains N-phenyl amide group\n'
               ' * SMILES: CCc1cccc(CC)c1N(COC\\C=C\\C)C(=O)CCl NAME: '
               'Butenachlor REASON: CORRECT Contains N-phenyl amide group\n'
               ' * SMILES: CC(C)N(C(=O)C(O)=O)c1ccccc1 NAME: propachlor OXA '
               'REASON: CORRECT Contains N-phenyl amide group\n'
               ' * SMILES: '
               'S([O-])(=O)(=O)CC(=O)N(C1=C(CC)C=CC=C1CC)COCCCC.[Na+] NAME: '
               'Butachlor ESA sodium salt REASON: CORRECT Contains N-phenyl '
               'amide group\n'
               'False positives: NONE\n'
               'False negatives: SMILES: '
               'C1=CC=C(C=C1)NC(=O)C[C@H]2C=C[C@H]([C@@H](O2)CO)NC(=O)CC3=CN=CC=C3 '
               'NAME: '
               '2-[(2R,3R,6S)-2-(hydroxymethyl)-3-[[1-oxo-2-(3-pyridinyl)ethyl]amino]-3,6-dihydro-2H-pyran-6-yl]-N-phenylacetamide '
               'REASON: MISSED No N-phenyl amide group found\n'
               ' * SMILES: '
               'CCN1C(=NN=C1SCC(=O)CC(=O)NC2=CC=CC=C2OC)COC3=CC=CC=C3C NAME: '
               '4-[[4-ethyl-5-[(2-methylphenoxy)methyl]-1,2,4-triazol-3-yl]thio]-N-(2-methoxyphenyl)-3-oxobutanamide '
               'REASON: MISSED No N-phenyl amide group found\n'
               ' * SMILES: '
               'CN(C)C=C(C=NOCC1=C(C=C(C=C1)Cl)Cl)C(=O)NC2=CC(=C(C=C2)Cl)Cl '
               'NAME: '
               'N-(3,4-dichlorophenyl)-2-[(2,4-dichlorophenyl)methoxyiminomethyl]-3-(dimethylamino)-2-propenamide '
               'REASON: MISSED No N-phenyl amide group found\n'
               ' * SMILES: '
               'CN1[C@@H]2CC[C@H](O[C@@H]2COC3=C(C1=O)C=C(C=C3)NC(=O)C4CCOCC4)CC(=O)NC5=CC(=CC=C5)OC '
               'NAME: '
               'N-[(2S,4aR,12aS)-2-[2-(3-methoxyanilino)-2-oxoethyl]-5-methyl-6-oxo-2,3,4,4a,12,12a-hexahydropyrano[2,3-c][1,5]benzoxazocin-8-yl]oxane-4-carboxamide '
               'REASON: MISSED No N-phenyl amide group found\n'
               ' * SMILES: N[C@@H](CCC(=O)Nc1ccccc1)C(O)=O NAME: '
               'N(5)-phenyl-L-glutamine REASON: MISSED No N-phenyl amide group '
               'found\n'
               ' * SMILES: COC1=CC=CC=C1NC(=O)CSC(=O)N NAME: carbamothioic '
               'acid S-[2-(2-methoxyanilino)-2-oxoethyl] ester REASON: MISSED '
               'No N-phenyl amide group found\n'
               ' * SMILES: '
               'CC1=CC(=NO1)C2=NC(=C(S2)C3=NN=C(N3C)SCC(=O)NC4=CC=C(C=C4)Cl)C '
               'NAME: '
               'N-(4-chlorophenyl)-2-[[4-methyl-5-[4-methyl-2-(5-methyl-3-isoxazolyl)-5-thiazolyl]-1,2,4-triazol-3-yl]thio]acetamide '
               'REASON: MISSED No N-phenyl amide group found\n'
               ' * SMILES: '
               'CC1(C(C1(C)C)C(=O)NC2=CC=C(C=C2)S(=O)(=O)C3=CC=C(C=C3)NC(=O)C4C(C4(C)C)(C)C)C '
               'NAME: '
               '2,2,3,3-tetramethyl-N-[4-[4-[[oxo-(2,2,3,3-tetramethylcyclopropyl)methyl]amino]phenyl]sulfonylphenyl]-1-cyclopropanecarboxamide '
               'REASON: MISSED No N-phenyl amide group found\n'
               ' * SMILES: C1CC1C(=O)NC2=CC=C(C=C2)NC3=CC=CC=C3 NAME: '
               'N-(4-anilinophenyl)cyclopropanecarboxamide REASON: MISSED No '
               'N-phenyl amide group found\n'
               ' * SMILES: CCCCC(=O)NN=C(C)C1=CC=C(C=C1)NC(=O)C(F)(F)F NAME: '
               'N-[1-[4-[(2,2,2-trifluoro-1-oxoethyl)amino]phenyl]ethylideneamino]pentanamide '
               'REASON: MISSED No N-phenyl amide group found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'N12CN3CN(CN(C1)CC3)CC2',
                                     'name': '1,3,6,8-tetraazatricyclo[4,4,1,1(3,8)]dodecane',
                                     'reason': 'No N-phenyl amide group found'},
                                 {   'smiles': 'CC1CC2C(CC=C1CCO)C(C)C(=O)O2',
                                     'name': '(8alpha,10beta,11beta)-3-hydroxy-4,15-dinor-1(5)-xanthen-12,8-olide',
                                     'reason': 'No N-phenyl amide group found'},
                                 {   'smiles': 'CCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](COP(O)(O)=O)[C@H](O)\\C=C\\CCCCCCCCCCCCC',
                                     'name': 'N-docosanoylsphingosine-1-phosphate',
                                     'reason': 'No N-phenyl amide group found'},
                                 {   'smiles': 'C1CNCCC1SC2=NC(=CC=C2)Cl',
                                     'name': '2-chloro-6-(4-piperidinylthio)pyridine',
                                     'reason': 'No N-phenyl amide group found'},
                                 {   'smiles': '[H][C@@]1(CC[C@@]2(C)[C@]3([H])CC[C@@]4([H])[C@@](C)([C@H](CC[C@@]44C[C@@]34CC[C@]12C)OC(C)=O)C(O)=O)[C@H](C)CCC(=O)C(C)C',
                                     'name': 'bonianic acid B, (rel)-',
                                     'reason': 'No N-phenyl amide group found'},
                                 {   'smiles': 'O=C1OC[C@@]2(C1=C(CC[C@]34[C@H]2C[C@H](CC[C@H]3C)C4(C)C)C)C',
                                     'name': 'Harzianelactone',
                                     'reason': 'No N-phenyl amide group found'},
                                 {   'smiles': 'O=C(NC1=NC=C(Cl)C=C1)CSC2=NN=C3C(NC4=C3C=CC=C4)=N2',
                                     'name': 'dCeMM2',
                                     'reason': 'No proper phenyl ring attached '
                                               'to amide nitrogen'},
                                 {   'smiles': 'N[C@@H](CS[C@H](\\C=C\\C=C\\C=C/C=C/CC(O)=O)[C@@H](O)CCCC(O)=O)C(O)=O',
                                     'name': '(13E)-16-carboxy-Delta(13)-17,18,19,20-tetranor-leukotriene '
                                             'E4',
                                     'reason': 'No N-phenyl amide group found'},
                                 {   'smiles': 'C[C@@H]1O[C@@H](OC[C@H]2OC(O)[C@H](NC(C)=O)[C@@H](O[C@@H]3O[C@@H](C)[C@@H](O)[C@@H](O)[C@@H]3O)[C@@H]2O[C@@H]2O[C@H](CO)[C@@H](O[C@@H]3O[C@H](CO[C@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]4O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4NC(C)=O)[C@@H](O)[C@H](O[C@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]4O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4NC(C)=O)[C@@H]3O)[C@H](O)[C@H]2NC(C)=O)[C@@H](O)[C@H](O)[C@@H]1O',
                                     'name': 'alpha-L-Fucp-(1->3)-[alpha-L-Fucp-(1->6)]-{beta-D-GlcpNAc-(1->2)-alpha-D-Manp-(1->3)-[beta-D-GlcpNAc-(1->2)-alpha-D-Manp-(1->6)]-beta-D-Manp-(1->4)-beta-GlcpNAc-(1->4)}-D-GlcpNAc',
                                     'reason': 'No N-phenyl amide group found'},
                                 {   'smiles': 'CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1ccc(N)nc1=O)OC(=O)CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CC',
                                     'name': 'CDP-1-stearoyl-2-(4Z,7Z,10Z,13Z,16Z,19Z)-docosahexaenoyl-sn-glycerol',
                                     'reason': 'No N-phenyl amide group '
                                               'found'}],
    'sample_false_negatives': [],
    'attempt': 3,
    'reasoning': 'Let me analyze the issues with the previous code:\n'
                 '\n'
                 '1. The main problem appears to be that the SMARTS pattern '
                 'was too restrictive:\n'
                 '- It required the amide carbon to NOT be bonded to N or O '
                 '`[!N;!O]`, which excluded valid anilides\n'
                 '- It required a specific substitution pattern on the '
                 'nitrogen `[NX3]([#1,#6])`\n'
                 '\n'
                 '2. Looking at the false negatives, we missed many valid '
                 'anilides where:\n'
                 '- The phenyl ring had substituents (like chloro, methoxy '
                 'groups)\n'
                 '- The nitrogen was part of a more complex structure\n'
                 '- The carbonyl carbon was connected to other groups\n'
                 '\n'
                 '3. The pattern should be simplified to just look for:\n'
                 '- A phenyl ring (can be substituted)\n'
                 '- Connected to a nitrogen \n'
                 '- That nitrogen is part of an amide group\n'
                 '\n'
                 "Here's my improved approach:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 249,
    'num_false_positives': 100,
    'num_true_negatives': 1025,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.7134670487106017,
    'recall': 1.0,
    'f1': 0.8327759197324415,
    'accuracy': 0.9272197962154294,
    'negative_predictive_value': 1.0}