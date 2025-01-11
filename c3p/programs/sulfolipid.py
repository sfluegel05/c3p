"""
Classifies: CHEBI:61384 sulfolipid
"""
"""
Classifies: sulfolipid
A compound containing a sulfate ester group attached to a lipid-containing sugar moiety.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sulfolipid(smiles: str):
    """
    Determines if a molecule is a sulfolipid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a sulfolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sulfate ester group (O-SO3H)
    sulfate_pattern = Chem.MolFromSmarts("[OX2][S](=[OX1])(=[OX1])[OX2H,OX1-]")
    if not mol.HasSubstructMatch(sulfate_pattern):
        return False, "No sulfate ester group found"

    # Common sulfolipid sugar patterns
    galactosyl_pattern = Chem.MolFromSmarts("[CH2]1O[CH]([CH]([CH]([CH]([CH]1O)O)O)O)O")
    trehalose_pattern = Chem.MolFromSmarts("O[CH]1[CH](O)[CH](O)[CH]([CH](O1)CO)O")
    
    has_galactosyl = mol.HasSubstructMatch(galactosyl_pattern)
    has_trehalose = mol.HasSubstructMatch(trehalose_pattern)
    
    if not (has_galactosyl or has_trehalose):
        return False, "Missing characteristic sugar moiety (galactosyl or trehalose)"

    # Look for lipid characteristics
    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:  # Increased threshold for complex sulfolipids
        return False, "Too few carbons for a sulfolipid structure"

    # Look for long carbon chains (at least 12 carbons)
    long_chain = Chem.MolFromSmarts("[CH2,CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]")
    if not mol.HasSubstructMatch(long_chain):
        return False, "No long carbon chains found"

    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:  # Increased threshold for complex sulfolipids
        return False, "Molecular weight too low for sulfolipid"

    # Look for common structural features
    ester_pattern = Chem.MolFromSmarts("[#6]C(=O)O[#6]")
    amide_pattern = Chem.MolFromSmarts("[#6]C(=O)N[#6]")
    sphingosine_pattern = Chem.MolFromSmarts("[CH2]O[CH]([CH]([CH]([CH]O)O)O)[CH]([CH](O)/C=C/[#6])N")
    
    features = []
    if mol.HasSubstructMatch(ester_pattern):
        features.append("ester bonds")
    if mol.HasSubstructMatch(amide_pattern):
        features.append("amide bonds")
    if mol.HasSubstructMatch(sphingosine_pattern):
        features.append("sphingosine backbone")
    if has_galactosyl:
        features.append("galactosyl moiety")
    if has_trehalose:
        features.append("trehalose moiety")

    # Build reason string
    reason = "Contains sulfate ester group with characteristic sugar-lipid structure"
    if features:
        reason += f" including {', '.join(features)}"

    return True, reason


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:61384',
                          'name': 'sulfolipid',
                          'definition': 'A compound containing a sulfonic acid '
                                        'residue joined by a carbon-sulfur '
                                        'bond to a lipid.',
                          'parents': ['CHEBI:18059', 'CHEBI:25704'],
                          'xrefs': ['Wikipedia:Sulfolipid'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: NONE\n'
               'False positives: SMILES: '
               'S(O)(=O)(=O)CC1OC(OCC(OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCCCCCCCCCC)C(O)C(O)C1O '
               'NAME: 1,2-Di-O-palmitoyl-3-O-(6-sulfoquinovopyranosyl)glycerol '
               'REASON: WRONGLY CLASSIFIED Contains sulfonic acid group with '
               'C-S bond and lipid characteristics including sugar moiety, '
               'ester bonds\n'
               ' * SMILES: '
               'S(O)(=O)(=O)CC(N)(C(OC(=O)CCCCCCCCCCCC(C)C)C(O)/C=C\\CCCCCCCCCCC(C)C)C(O)=O '
               'NAME: '
               "2-carboxy-2-amino-3-O-(13'-methyltetradecanoyl)-4-hydroxy-17-methyloctadec-5-ene-1-sulfonic "
               'acid REASON: WRONGLY CLASSIFIED Contains sulfonic acid group '
               'with C-S bond and lipid characteristics including ester bonds\n'
               ' * SMILES: C(C(NCCS(O)(=O)=O)=O)CCCCCCCCCCCC/C=C\\CCCCCCCC '
               'NAME: N-nervonoyltaurine REASON: WRONGLY CLASSIFIED Contains '
               'sulfonic acid group with C-S bond and lipid characteristics '
               'including amide bonds\n'
               ' * SMILES: C(C(NCCS(O)(=O)=O)=O)CCCCCC/C=C\\CCCCCCCC NAME: '
               'N-oleoyltaurine REASON: WRONGLY CLASSIFIED Contains sulfonic '
               'acid group with C-S bond and lipid characteristics including '
               'amide bonds\n'
               ' * SMILES: CC(CCCCCCCCCCCC([C@@H]([NH3+])CS(=O)(=O)[O-])O)C '
               'NAME: 2-amino-3-hydroxy-15-methylhexadecane-1-sulfonate '
               'REASON: WRONGLY CLASSIFIED Contains sulfonic acid group with '
               'C-S bond and lipid characteristics\n'
               ' * SMILES: CCCCCCCCCCCCCCCC([C@@H]([NH3+])CS(=O)(=O)[O-])=O '
               'NAME: 2-amino-3-oxo-octadecane-1-sulfonate REASON: WRONGLY '
               'CLASSIFIED Contains sulfonic acid group with C-S bond and '
               'lipid characteristics\n'
               ' * SMILES: '
               '[C@H]1(O)C(OCC(OC(CCCCCCCCCCCCC)=O)COC(=O)CCCCCCCCCCCCC)O[C@H](CS(O)(=O)=O)[C@H]([C@@H]1O)O '
               'NAME: 1,2-ditetradecanoyl-3-(6-sulfoquinovopyranosyl)glycerol '
               'REASON: WRONGLY CLASSIFIED Contains sulfonic acid group with '
               'C-S bond and lipid characteristics including sugar moiety, '
               'ester bonds\n'
               ' * SMILES: '
               'S(O)(=O)(=O)C[C@H]1O[C@H](OC[C@H](OC(=O)CCCCCCCCCCC/C=C\\CC)CO)C(O)C(O)[C@@H]1O '
               'NAME: SQMG(0:0/16:1(13Z)) REASON: WRONGLY CLASSIFIED Contains '
               'sulfonic acid group with C-S bond and lipid characteristics '
               'including sugar moiety, ester bonds\n'
               ' * SMILES: C(C(NCCS(O)(=O)=O)=O)CCCCCCCCCCCCCCCCCC NAME: '
               'N-icosanoyltaurine REASON: WRONGLY CLASSIFIED Contains '
               'sulfonic acid group with C-S bond and lipid characteristics '
               'including amide bonds\n'
               ' * SMILES: CCCCCCCCCCCCCC[N+](C)(C)CCCS([O-])(=O)=O NAME: '
               'Zwittergent 3-14 REASON: WRONGLY CLASSIFIED Contains sulfonic '
               'acid group with C-S bond and lipid characteristics\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCC(=O)OCC(COP([O-])(=O)OCCNS(=O)(=O)c1ccc(c(c1)S([O-])(=O)=O)C1=c2cc3CCC[N+]4=c3c(CCC4)c2Oc2c3CCCN4CCCc(cc12)c34)OC(=O)CCCCCCCCCCCCCCC '
               'NAME: texas red DHPE(1-) REASON: WRONGLY CLASSIFIED Contains '
               'sulfonic acid group with C-S bond and lipid characteristics '
               'including ester bonds\n'
               ' * SMILES: C(C(NCCS([O-])(=O)=O)=O)CCCCCCCCCCCCCCCCCCCCCC '
               'NAME: N-tetracosanoyltaurine(1-) REASON: WRONGLY CLASSIFIED '
               'Contains sulfonic acid group with C-S bond and lipid '
               'characteristics including amide bonds\n'
               ' * SMILES: C(C(NCCS(O)(=O)=O)=O)CCCCCCCCCCCCCCCC NAME: '
               'N-stearoyltaurine REASON: WRONGLY CLASSIFIED Contains sulfonic '
               'acid group with C-S bond and lipid characteristics including '
               'amide bonds\n'
               ' * SMILES: '
               'S(=O)(=O)(O)C[C@H](NC(=O)[C@H](O)CCCCCCCCCCC(C)C)[C@H](O)C[C@H](O)CCCCCCCCCC(C)C '
               'NAME: RIF-1 REASON: WRONGLY CLASSIFIED Contains sulfonic acid '
               'group with C-S bond and lipid characteristics including amide '
               'bonds\n'
               ' * SMILES: C(C(NCCS([O-])(=O)=O)=O)CCCCCCCCCCCCCCCCCC NAME: '
               'N-icosanoyltaurine(1-) REASON: WRONGLY CLASSIFIED Contains '
               'sulfonic acid group with C-S bond and lipid characteristics '
               'including amide bonds\n'
               ' * SMILES: '
               'S(=O)(=O)(O)C[C@H](NC(=O)C[C@H](O)CCCCCCCCCCCC(C)C)[C@H](O)CCCCCCCCCCCC(C)C '
               'NAME: Sulfobacin A REASON: WRONGLY CLASSIFIED Contains '
               'sulfonic acid group with C-S bond and lipid characteristics '
               'including amide bonds\n'
               ' * SMILES: S(O)(=O)(=O)CCNC(=O)CCCCCCCCCCCCCCCCCCCCCCCC NAME: '
               'n-pentacosanoyl taurine REASON: WRONGLY CLASSIFIED Contains '
               'sulfonic acid group with C-S bond and lipid characteristics '
               'including amide bonds\n'
               ' * SMILES: '
               'S(O)(=O)(=O)C[C@H](NC(=O)CCCCCCCCCCCC(C)C)[C@H](O)CCCCCCCCCCCC(C)C '
               'NAME: N-(13-methyl-tetradecanoyl)-capnine REASON: WRONGLY '
               'CLASSIFIED Contains sulfonic acid group with C-S bond and '
               'lipid characteristics including amide bonds\n'
               ' * SMILES: '
               'S(O)(=O)(=O)CCOC(=O)/C=C(\\CCCCC[C@H](OC(=O)CCCCCCCCCCCCC/C=C\\CCCCCC)C)/C '
               'NAME: '
               '2-(9R-(15Z-docosenoyloxy)-3-methyl-2Z-decenoyloxy)-ethanesulfonic '
               'acid REASON: WRONGLY CLASSIFIED Contains sulfonic acid group '
               'with C-S bond and lipid characteristics including ester bonds\n'
               ' * SMILES: '
               'S(=O)(=O)(O)CC1OC(OCC(OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCC/C=C\\CCCCCCCC)C(OC(=O)CCCCCCCCCCCCCCC)C(C1O)OC(=O)CCCCCCCCCCCCCCC '
               'NAME: '
               '(3-hydroxy-6-(3-(oleoyloxy)-2-(palmitoyloxy)propoxy)-4,5-bis(palmitoyloxy)tetrahydro-2H-pyran-2-yl)methanesulfonic '
               'acid REASON: WRONGLY CLASSIFIED Contains sulfonic acid group '
               'with C-S bond and lipid characteristics including sugar '
               'moiety, ester bonds\n'
               ' * SMILES: '
               'S(O)(=O)(=O)C[C@H]1O[C@H](OC[C@H](O)COC(=O)CCCCCCCCC/C=C\\CCCC)C(O)C(O)[C@@H]1O '
               'NAME: SQMG(16:1(11Z)/0:0) REASON: WRONGLY CLASSIFIED Contains '
               'sulfonic acid group with C-S bond and lipid characteristics '
               'including sugar moiety, ester bonds\n'
               ' * SMILES: '
               'S(=O)(=O)(O)C[C@@H](NC(=O)C[C@H](O)CCCCCCCCCCCC(C)C)[C@H](O)/C=C/CCCCCCCCCC(C)C '
               'NAME: Flavocristamide A REASON: WRONGLY CLASSIFIED Contains '
               'sulfonic acid group with C-S bond and lipid characteristics '
               'including amide bonds\n'
               ' * SMILES: '
               'S(O)(=O)(=O)C[C@H](NC(=O)C\\C=C\\CCCCCCCCCCCC)[C@H](O)\\C=C\\CCCCCCCCCCCCC '
               'NAME: N-(3E-hexadecenoyl)-deoxysphing-4-enine-1-sulfonate '
               'REASON: WRONGLY CLASSIFIED Contains sulfonic acid group with '
               'C-S bond and lipid characteristics including amide bonds\n'
               ' * SMILES: C(C(NCCS(O)(=O)=O)=O)CCCCCCCCCCCCCCCCCCCC NAME: '
               'N-docosanoyltaurine REASON: WRONGLY CLASSIFIED Contains '
               'sulfonic acid group with C-S bond and lipid characteristics '
               'including amide bonds\n'
               ' * SMILES: C(C(NCCS([O-])(=O)=O)=O)CCCCCCCCCCCCCCCCCCCCC NAME: '
               'N-tricosanoyltaurine(1-) REASON: WRONGLY CLASSIFIED Contains '
               'sulfonic acid group with C-S bond and lipid characteristics '
               'including amide bonds\n'
               'False negatives: SMILES: '
               'S(OC1[C@@H](O)[C@H](O[C@@H](OC[C@H](NC(=O)CCCCCCCCCCCCCCC/C=C\\CCCCCCCC)[C@H](O)/C=C/CCCCCCCCCCCCC)C1O)CO)(O)(=O)=O '
               'NAME: '
               '(17Z)-N-[(2S,3R,4E)-1-{[(2R,5S,6R)-3,5-dihydroxy-6-(hydroxymethyl)-4-(sulfooxy)oxan-2-yl]oxy}-3-hydroxyoctadec-4-en-2-yl]hexacos-17-enimidic '
               'acid REASON: MISSED No carbon-sulfur bond to sulfonic acid '
               'group found\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCC[C@H](C)\\C=C(/C)C(=O)O[C@H]1[C@H](O)[C@@H](CO)O[C@H](O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2OS(O)(=O)=O)[C@@H]1OC(=O)CCCCCCCCCCCCCCC '
               'NAME: '
               "3-O-[(2E,4S)-2,4-dimethyldocos-2-enoyl]-2-O-palmitoyl-2'-O-sulfo-alpha,alpha-trehalose "
               'REASON: MISSED No carbon-sulfur bond to sulfonic acid group '
               'found\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](OS(O)(=O)=O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC '
               'NAME: 1-(3-O-sulfo-beta-D-galactosyl)-N-palmitoylsphingosine '
               'REASON: MISSED No carbon-sulfur bond to sulfonic acid group '
               'found\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCC[C@H](C)C[C@H](C)\\C=C(/C)C(=O)O[C@H]1[C@H](O)[C@@H](CO)O[C@H](O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2OS(O)(=O)=O)[C@@H]1OC(=O)CCCCCCCCCCCCCCC '
               'NAME: '
               "2-O-palmitoyl-3-O-[(2E,4S,6S)-2,4,6-trimethyltetracos-2-enoyl]-2'-O-sulfo-alpha,alpha-trehalose "
               'REASON: MISSED No carbon-sulfur bond to sulfonic acid group '
               'found\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCC(=O)O[C@H]1[C@H](O[C@H](CO)[C@@H](O)[C@@H]1OC(=O)CCCCCCCCCCCCCCCO)O[C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1OS(O)(=O)=O '
               'NAME: '
               "3-O-(16-hydroxyhexadecanoyl)-2-O-palmitoyl-2'-O-sulfo-alpha,alpha-trehalose "
               'REASON: MISSED No carbon-sulfur bond to sulfonic acid group '
               'found\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCC(C)C(=O)O[C@H]1[C@H](O)[C@@H](CO)O[C@H](O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2OS(O)(=O)=O)[C@@H]1OC(=O)CCCCCCCCCCCCCCC '
               'NAME: '
               "3-O-(2-methyloctadecanoyl)-2-O-palmitoyl-2'-O-sulfo-alpha,alpha-trehalose "
               'REASON: MISSED No carbon-sulfur bond to sulfonic acid group '
               'found\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](OS(O)(=O)=O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC '
               'NAME: 1-(3-O-sulfo-beta-D-galactosyl)-N-docosanoylsphingosine '
               'REASON: MISSED No carbon-sulfur bond to sulfonic acid group '
               'found\n'
               ' * SMILES: '
               '[C@H]1([C@@H]([C@H]([C@@H]([C@H](O1)CO)O)O)OS(=O)(=O)O)O[C@@H]2[C@@H]([C@H]([C@@H]([C@H](O2)CO)O)OC(/C(=C/[C@H](C[C@H](C[C@H](CCCCCCCCCCCCCCCCCCCCCC)C)C)C)/C)=O)OC(CCCCCCCCCCCCCCC)=O '
               'NAME: '
               "2-O-palmitoyl-2'-O-sulfo-3-O-[(2E,4S,6S,8S)-2,4,6,8-tetramethyltriacont-2-enoyl]-alpha,alpha-trehalose "
               'REASON: MISSED No carbon-sulfur bond to sulfonic acid group '
               'found\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCC[C@H](C)C[C@H](C)C(=O)O[C@H]1[C@H](O)[C@@H](CO)O[C@H](O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2OS(O)(=O)=O)[C@@H]1OC(=O)CCCCCCCCCCCCCCC '
               'NAME: '
               "3-O-[(2S,4S)-2,4-dimethyldocosanoyl]-2-O-palmitoyl-2'-O-sulfo-alpha,alpha-trehalose "
               'REASON: MISSED No carbon-sulfur bond to sulfonic acid group '
               'found\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCCCCCCCC(=O)O[C@H]1[C@H](O[C@H](CO)[C@@H](O)[C@@H]1OC(=O)C(\\C)=C\\[C@@H](C)C[C@@H](C)C[C@@H](C)CCCCCCCCCCCCCCCC)O[C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1OS(O)(=O)=O '
               'NAME: '
               "2'-O-sulfo-2-O-tetracosanoyl-3-O-[(2E,4S,6S,8S)-2,4,6,8-tetramethyltetracos-2-enoyl]-alpha,alpha-trehalose "
               'REASON: MISSED No carbon-sulfur bond to sulfonic acid group '
               'found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C1C=C(N=C2N1NN=N2)C3=CC=CC=C3',
                                     'name': '5-phenyl-1,7-dihydrotetrazolo[1,5-a]pyrimidine',
                                     'reason': 'No sulfate ester group found'},
                                 {   'smiles': 'O(C(=O)CCCCCCCCCCC/C=C\\CCCCCCCC)[C@H](COC(=O)CCCCCCCCC/C=C\\CCCCCC)COC(=O)CCCCCCC/C=C\\CCCC',
                                     'name': 'TG(14:1(9Z)/22:1(13Z)/18:1(11Z))',
                                     'reason': 'No sulfate ester group found'},
                                 {   'smiles': 'CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@@H](O[C@@H]3O[C@H](CO[C@H]4O[C@H](CO[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H](O)[C@H](O[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H]4O)[C@@H](O)[C@H](O[C@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]4O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4NC(C)=O)[C@@H]3O)[C@H](O)[C@H]2NC(C)=O)[C@@H]1O',
                                     'name': 'beta-D-GlcpNAc-(1->2)-alpha-D-Manp-(1->3)-{alpha-D-Manp-(1->3)-[alpha-D-Manp-(1->6)]-alpha-D-Manp-(1->6)}-beta-D-Manp-(1->4)-beta-GlcpNAc-(1->4)-beta-D-GlcpNAc',
                                     'reason': 'No sulfate ester group found'},
                                 {   'smiles': '[C@@]12([C@@]([C@]([C@@H](CC1)C)(CCC(CC)C)C)(CCC[C@@H]2C)[H])C',
                                     'name': 'clerodane',
                                     'reason': 'No sulfate ester group found'},
                                 {   'smiles': 'S(O[C@@H]1[C@H](O)[C@@H](O)[C@H](O[C@@H]([C@@H](O)[C@H](O)CO[C@]2(O[C@H]([C@H](NC(=O)C)[C@@H](O)C2)[C@H](O)[C@H](O)CO)C(O)=O)[C@@H](NC(=O)C)CO)O[C@@H]1CO)(O)(=O)=O',
                                     'name': '(2R,4S,5R,6R)-5-Acetamido-2-[(2R,3S,4R,5S)-5-acetamido-4-[(2R,3R,4R,5R,6R)-3,4-dihydroxy-6-(hydroxymethyl)-5-sulfooxyoxan-2-yl]oxy-2,3,6-trihydroxyhexoxy]-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'No long carbon chains found'},
                                 {   'smiles': 'O=C(NC(CC(C)C)C(=O)N[C@@H](CCCN=C(N)N)C=O)C(NC(=O)CC)CC(C)C',
                                     'name': 'Leupeptin Pr-LL',
                                     'reason': 'No sulfate ester group found'},
                                 {   'smiles': 'COC1=CC=C(C=C1)CCNC(=O)C(=C2C3=CC=CC=C3C(=N2)NC(=O)C4=CC=CS4)C#N',
                                     'name': 'N-[3-[1-cyano-2-[2-(4-methoxyphenyl)ethylamino]-2-oxoethylidene]-1-isoindolyl]-2-thiophenecarboxamide',
                                     'reason': 'No sulfate ester group found'},
                                 {   'smiles': 'O=C1O[C@@H](CC=2C1=C(OC)C(O)=CC2)CCC[C@@H](O)C',
                                     'name': 'Penicimarin C',
                                     'reason': 'No sulfate ester group found'},
                                 {   'smiles': 'CO[C@@H]1[C@@H]2[C@@H](C[C@@H]3[C@@](O2)([C@H]([C@@H]([C@H](O3)C(=O)OC)O)O)O)OC1N4C=NC5=C4N=CN=C5N',
                                     'name': 'LSM-4497',
                                     'reason': 'No sulfate ester group found'},
                                 {   'smiles': 'C[C@@H]1[C@H](O)CCC2=CC[C@H](C[C@]12C)C(C)=C',
                                     'name': '1-deoxycapsidiol',
                                     'reason': 'No sulfate ester group found'}],
    'sample_false_negatives': [   {   'smiles': 'C(=C/CCCCCCCCCCCCC)\\[C@@H](O)[C@@H](N)CO[C@H]1[C@@H]([C@H]([C@H]([C@H](O1)CO)O)OS(O)(=O)=O)O',
                                      'name': '1-(3-O-sulfo-beta-D-galactosyl)sphingosine',
                                      'reason': 'Molecular weight too low for '
                                                'sulfolipid'},
                                  {   'smiles': 'CCCCCCCCCCCCC\\C=C\\[C@@H](O)[C@@H](N)CO[C@@H]1O[C@H](COS(O)(=O)=O)[C@H](O)[C@H](O)[C@H]1O',
                                      'name': 'psychosine sulfate',
                                      'reason': 'Molecular weight too low for '
                                                'sulfolipid'}],
    'attempt': 2,
    'reasoning': 'Looking at the outcomes, I can see several key issues with '
                 'the previous implementation:\n'
                 '\n'
                 "1. The main problem is that we're misinterpreting sulfate "
                 'esters as sulfonic acids. The true positive examples all '
                 'contain sulfate esters (O-SO3H groups) rather than sulfonic '
                 'acids (C-SO3H groups). The definition appears to be '
                 'incorrect or misleading.\n'
                 '\n'
                 "2. We're getting many false positives with actual sulfonic "
                 'acids (taurine derivatives, etc.) because we were looking '
                 'for C-S bonds.\n'
                 '\n'
                 "3. We're missing all the true sulfolipids because they "
                 'contain sulfate esters (O-SO3H) attached to sugar or lipid '
                 'moieties.\n'
                 '\n'
                 "Here's how we should improve the implementation:\n"
                 '\n'
                 '1. Change the primary search pattern to look for sulfate '
                 'esters (O-SO3H) instead of sulfonic acids (C-SO3H)\n'
                 '2. Add more specific patterns to identify the sugar-lipid '
                 'components that are characteristic of these molecules\n'
                 '3. Strengthen the lipid requirements to avoid false '
                 'positives\n'
                 '4. Add specific patterns for common sulfolipid structural '
                 'features like galactosyl and trehalose\n'
                 '\n'
                 "Here's the corrected implementation:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 32,
    'num_false_positives': 26,
    'num_true_negatives': 142240,
    'num_false_negatives': 2,
    'num_negatives': None,
    'precision': 0.5517241379310345,
    'recall': 0.9411764705882353,
    'f1': 0.6956521739130435,
    'accuracy': 0.999803232607168,
    'negative_predictive_value': 0.9999859394552945}