"""
Classifies: CHEBI:52221 isothiocyanate
"""
"""
Classifies: isothiocyanate compounds (R-N=C=S)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_isothiocyanate(smiles: str):
    """
    Determines if a molecule is an isothiocyanate based on its SMILES string.
    Isothiocyanates have the general formula R-N=C=S where R is an organic group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isothiocyanate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for N=C=S group
    isothiocyanate_pattern = Chem.MolFromSmarts("[NX2]=[CX2]=[SX1]")
    matches = mol.GetSubstructMatches(isothiocyanate_pattern)
    
    if not matches:
        return False, "No isothiocyanate group found"
    
    # Get basic molecular properties
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    num_atoms = mol.GetNumAtoms()
    
    valid_groups = 0
    for match in matches:
        n_atom, c_atom, s_atom = [mol.GetAtomWithIdx(idx) for idx in match]
        
        # Validate nitrogen
        if n_atom.GetDegree() != 2 or n_atom.GetTotalValence() != 3:
            continue
            
        # Validate carbon
        if c_atom.GetDegree() != 2 or c_atom.GetTotalValence() != 4:
            continue
            
        # Validate sulfur
        if s_atom.GetDegree() != 1 or s_atom.GetTotalValence() != 2:
            continue
            
        # Get R group (atom connected to nitrogen that's not the central carbon)
        r_atoms = [atom for atom in n_atom.GetNeighbors() if atom.GetIdx() != c_atom.GetIdx()]
        if not r_atoms:
            continue
        r_atom = r_atoms[0]
        
        # Validate R group
        if r_atom.GetAtomicNum() not in [6, 7, 8]:  # Must be C, N, or O
            continue
            
        # Check bond types
        n_c_bond = mol.GetBondBetweenAtoms(n_atom.GetIdx(), c_atom.GetIdx())
        c_s_bond = mol.GetBondBetweenAtoms(c_atom.GetIdx(), s_atom.GetIdx())
        
        if not (n_c_bond.GetBondType() == Chem.BondType.DOUBLE and 
                c_s_bond.GetBondType() == Chem.BondType.DOUBLE):
            continue
        
        # Check if the isothiocyanate group is a significant feature
        # by comparing its size to the overall molecule
        ncs_size = 3  # N=C=S group size
        if ncs_size / num_atoms < 0.1:  # If group is less than 10% of molecule
            if mol_weight > 500:  # And molecule is large
                continue  # Likely a complex molecule where -NCS is not the main feature
        
        valid_groups += 1
    
    if valid_groups == 0:
        return False, "No valid isothiocyanate groups found"
        
    # Additional checks for problematic features
    problematic_patterns = [
        Chem.MolFromSmarts("[N+]"), # Quaternary nitrogen
        Chem.MolFromSmarts("[N-]"), # Negatively charged nitrogen
        Chem.MolFromSmarts("[S+]"), # Positively charged sulfur
        Chem.MolFromSmarts("[S-]"), # Negatively charged sulfur
    ]
    
    for pattern in problematic_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            return False, "Contains charged species that indicate different chemical class"
    
    # Check if molecule is too complex
    if mol_weight > 800:
        return False, "Molecule too complex to be classified primarily as isothiocyanate"
    
    if valid_groups == 1:
        return True, "Contains one isothiocyanate group (R-N=C=S)"
    else:
        return True, f"Contains {valid_groups} isothiocyanate groups"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:52221',
                          'name': 'isothiocyanate',
                          'definition': 'An organosulfur compound with the '
                                        'general formula R-N=C=S.',
                          'parents': ['CHEBI:33261'],
                          'xrefs': [   'MetaCyc:ISOTHIOCYANATES',
                                       'PMID:21783213',
                                       'Wikipedia:Isothiocyanate'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.62 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: C=CCN=C=S NAME: allyl isothiocyanate '
               'REASON: CORRECT Contains one isothiocyanate group (R-N=C=S)\n'
               ' * SMILES: N(CCCS(C)=O)=C=S NAME: iberin REASON: CORRECT '
               'Contains one isothiocyanate group (R-N=C=S)\n'
               ' * SMILES: CCN=C=S NAME: ethyl isothiocyanate REASON: CORRECT '
               'Contains one isothiocyanate group (R-N=C=S)\n'
               ' * SMILES: Cc1ccc(cc1N=C=S)N=C=S NAME: toluene '
               '2,4-diisothiocyanate REASON: CORRECT Contains 2 isothiocyanate '
               'groups\n'
               ' * SMILES: C=CCCCN=C=S NAME: 5-isothiocyanato-1-pentene '
               'REASON: CORRECT Contains one isothiocyanate group (R-N=C=S)\n'
               ' * SMILES: N(=C=S)CCCC[S@@](=O)C NAME: (S)-sulforaphane '
               'REASON: CORRECT Contains one isothiocyanate group (R-N=C=S)\n'
               ' * SMILES: S(CCCCN=C=S)C NAME: Erucin REASON: CORRECT Contains '
               'one isothiocyanate group (R-N=C=S)\n'
               ' * SMILES: S=C=Nc1cccc2ccccc12 NAME: 1-naphthyl isothiocyanate '
               'REASON: CORRECT Contains one isothiocyanate group (R-N=C=S)\n'
               ' * SMILES: C1CC1N=C=S NAME: Isothiocyanatocyclopropane REASON: '
               'CORRECT Contains one isothiocyanate group (R-N=C=S)\n'
               ' * SMILES: N(=C=S)CCCC[S@](=O)C NAME: (R)-sulforaphane REASON: '
               'CORRECT Contains one isothiocyanate group (R-N=C=S)\n'
               ' * SMILES: OCCCN=C=S NAME: 3-hydroxypropyl isothiocyanate '
               'REASON: CORRECT Contains one isothiocyanate group (R-N=C=S)\n'
               ' * SMILES: S=C=Nc1ccccc1 NAME: phenyl isothiocyanate REASON: '
               'CORRECT Contains one isothiocyanate group (R-N=C=S)\n'
               ' * SMILES: COC1=CC=CC2=C1C(=CN2)CN=C=S NAME: '
               '(4-methoxyindol-3-ylmethyl)isothiocyanate REASON: CORRECT '
               'Contains one isothiocyanate group (R-N=C=S)\n'
               ' * SMILES: S=C=NCCc1ccccc1 NAME: phenethyl isothiocyanate '
               'REASON: CORRECT Contains one isothiocyanate group (R-N=C=S)\n'
               ' * SMILES: C(CCCCS(C)=O)CCCN=C=S NAME: 8-(methylsulfinyl)octyl '
               'isothiocyanate REASON: CORRECT Contains one isothiocyanate '
               'group (R-N=C=S)\n'
               ' * SMILES: '
               '[Cl-].CCN(CC)c1ccc2c(-c3ccc(cc3C(O)=O)N=C=S)c3ccc(cc3oc2c1)=[N+](CC)CC '
               'NAME: rhodamine B 5-isothiocyanate REASON: CORRECT Contains '
               'one isothiocyanate group (R-N=C=S)\n'
               ' * SMILES: CN=C=S NAME: methyl isothiocyanate REASON: CORRECT '
               'Contains one isothiocyanate group (R-N=C=S)\n'
               ' * SMILES: '
               'CCN(CC)c1ccc2c(-c3cc(ccc3C(O)=O)N=C=S)c3ccc(cc3oc2c1)=[N+](CC)CC '
               'NAME: rhodamine B 6-isothiocyanate REASON: CORRECT Contains '
               'one isothiocyanate group (R-N=C=S)\n'
               ' * SMILES: S(CCCCCN=C=S)C NAME: Berteroin REASON: CORRECT '
               'Contains one isothiocyanate group (R-N=C=S)\n'
               ' * SMILES: S=C=NC1CCCCC1 NAME: isothiocyanato-Cyclohexane '
               'REASON: CORRECT Contains one isothiocyanate group (R-N=C=S)\n'
               ' * SMILES: N(=C=S)CCCCS(=O)C NAME: sulforaphane REASON: '
               'CORRECT Contains one isothiocyanate group (R-N=C=S)\n'
               ' * SMILES: Cc1c(cccc1N=C=S)N=C=S NAME: toluene '
               '2,6-diisothiocyanate REASON: CORRECT Contains 2 isothiocyanate '
               'groups\n'
               ' * SMILES: S(CCCCCCCCN=C=S)C NAME: '
               '8-(methylthio)octylisothiocyanate REASON: CORRECT Contains one '
               'isothiocyanate group (R-N=C=S)\n'
               ' * SMILES: S(SCCCCN=C=S)CCCCN=C=S NAME: '
               'Bis(4-isothiocyanatobutyl) disulfide REASON: CORRECT Contains '
               '2 isothiocyanate groups\n'
               ' * SMILES: S(CCCCCCCN=C=S)C NAME: '
               '1-Isothiocyanato-7-(methylthio)heptane REASON: CORRECT '
               'Contains one isothiocyanate group (R-N=C=S)\n'
               'False positives: SMILES: S=C=NCC1=CC=C(O)C=C1 NAME: '
               '4-Hydroxybenzyl isothiocyanate REASON: WRONGLY CLASSIFIED '
               'Contains one isothiocyanate group (R-N=C=S)\n'
               ' * SMILES: S=C=NCCCCC1=CC=CC=C1 NAME: '
               '1-Isothiocyanato-4-phenylbutane REASON: WRONGLY CLASSIFIED '
               'Contains one isothiocyanate group (R-N=C=S)\n'
               ' * SMILES: '
               'OS(=O)(=O)c1cc(ccc1\\C=C\\c1ccc(cc1S(O)(=O)=O)N=C=S)N=C=S '
               "NAME: 4,4'-diisothiocyano-trans-stilbene-2,2'-disulfonic acid "
               'REASON: WRONGLY CLASSIFIED Contains 2 isothiocyanate groups\n'
               ' * SMILES: '
               'S=C=N[C@H]1[C@@](C=C)(CC[C@H]([C@@H]1C=2C3=C(C=CC=C3)NC2)C(=C)C)C '
               'NAME: Hapalindole Q REASON: WRONGLY CLASSIFIED Contains one '
               'isothiocyanate group (R-N=C=S)\n'
               ' * SMILES: Oc1ccc2c(Oc3cc(O)ccc3C22OC(=O)c3ccc(cc23)N=C=S)c1 '
               'NAME: fluorescein 6-isothiocyanate REASON: WRONGLY CLASSIFIED '
               'Contains one isothiocyanate group (R-N=C=S)\n'
               ' * SMILES: '
               'S=C=N[C@@H]1C2=C3C(=O)NC=4C3=C(C=CC4)C([C@H]2CC[C@@]1(C=C)C)(C)C '
               'NAME: Anhydrohapaloxindole M REASON: WRONGLY CLASSIFIED '
               'Contains one isothiocyanate group (R-N=C=S)\n'
               ' * SMILES: '
               '[O-]Cl(=O)(=O)=O.CN(C)c1ccc(cc1)C(c1ccc(cc1)N=C=S)=C1C=CC(C=C1)=[N+](C)C '
               'NAME: malachite green isothiocyanate REASON: WRONGLY '
               'CLASSIFIED Contains one isothiocyanate group (R-N=C=S)\n'
               ' * SMILES: '
               'Cl[C@@H]1C(O[C@]([C@@H]2[C@H]3[C@@H](NC=O)[C@](O)(CC[C@@H]3[C@@](N=C=S)(C)CC2)C)(C)CC1)(C)C '
               'NAME: Kalihinol J REASON: WRONGLY CLASSIFIED Contains one '
               'isothiocyanate group (R-N=C=S)\n'
               ' * SMILES: '
               'Cl[C@H]1[C@@](C=C)([C@@]2(N=C=S)C3=C4C(NC([C@@H]4C([C@@H](C1)C2=O)(C)C)=O)=CC=C3)C '
               'NAME: Welwitindolinone B isothiocyanate REASON: WRONGLY '
               'CLASSIFIED Contains one isothiocyanate group (R-N=C=S)\n'
               ' * SMILES: S=C=NC1=CC=C(C=C1)N=C=S NAME: bitoscanate REASON: '
               'WRONGLY CLASSIFIED Contains 2 isothiocyanate groups\n'
               ' * SMILES: '
               'S=C=N/C(/C(=O)O[C@@H]1[C@@H](O[C@H]2O[C@H]([C@](O)([C@@H](OC(=O)C(C)C)C)[C@@H](C2)O)C)[C@H](O)[C@H]([C@@]3(O)C(O)=C(C(=N)C(C3)=O)C(=O)O)O[C@H]1COC(=O)C)=C/C '
               'NAME: O-demethylpaulomycin B REASON: WRONGLY CLASSIFIED '
               'Contains one isothiocyanate group (R-N=C=S)\n'
               ' * SMILES: '
               'S=C=N/C(/C(=O)O[C@H]1[C@@H](O[C@@H]2O[C@@H]([C@@](O)([C@@H](O)C)[C@@H](C2)OC)C)[C@@H](O)[C@H]([C@@]3(O)C(O)=C(C(=N)C(C3)=O)C(=O)O)O[C@H]1COC(=O)C)=C\\C '
               'NAME: Paulomycin F REASON: WRONGLY CLASSIFIED Contains one '
               'isothiocyanate group (R-N=C=S)\n'
               ' * SMILES: '
               'CC(=O)Nc1ccc(\\C=C\\c2ccc(cc2S(O)(=O)=O)N=C=S)c(c1)S(O)(=O)=O '
               "NAME: 4-Acetamido-4'-isothiocyanostilbene-2,2'-disulphonic "
               'acid REASON: WRONGLY CLASSIFIED Contains one isothiocyanate '
               'group (R-N=C=S)\n'
               ' * SMILES: '
               'OS(=O)(=O)c1cc(ccc1\\C=C/c1ccc(cc1S(O)(=O)=O)N=C=S)N=C=S NAME: '
               "4,4'-diisothiocyano-cis-stilbene-2,2'-disulfonic acid REASON: "
               'WRONGLY CLASSIFIED Contains 2 isothiocyanate groups\n'
               ' * SMILES: [O-]C(=O)c1ccc(cc1)N=C=S NAME: '
               '4-isothiocyanatobenzoate REASON: WRONGLY CLASSIFIED Contains '
               'one isothiocyanate group (R-N=C=S)\n'
               ' * SMILES: '
               'Cl[C@H]1[C@](C=C)([C@H](N=C=S)[C@H]2C=3C=4C(=CC=CC4C([C@H]2C1)(C)C)NC3)C '
               'NAME: Hapalindole B REASON: WRONGLY CLASSIFIED Contains one '
               'isothiocyanate group (R-N=C=S)\n'
               ' * SMILES: S(=O)(CCCCCCN=C=S)C NAME: '
               '1-isothiocyanato-6-(methylsulfinyl)hexane REASON: WRONGLY '
               'CLASSIFIED Contains one isothiocyanate group (R-N=C=S)\n'
               ' * SMILES: '
               'S=C=N/C(/C(=O)O[C@H]1[C@@H](O[C@@H]2O[C@H]([C@](O)([C@H](OC(=O)[C@H](CC)C)C)[C@H](C2)OC)C)[C@@H](O)[C@H]([C@]3(O)C(O)=C(C(=N)C(C3)=O)C(=O)O)O[C@@H]1COC(=O)C)=C/C '
               'NAME: Paulomycin A REASON: WRONGLY CLASSIFIED Contains one '
               'isothiocyanate group (R-N=C=S)\n'
               ' * SMILES: S(=O)(=O)(CCCCCCN=C=S)C NAME: '
               '1-Isothiocyanato-6-(methylsulfonyl)-hexane REASON: WRONGLY '
               'CLASSIFIED Contains one isothiocyanate group (R-N=C=S)\n'
               ' * SMILES: '
               'S=C=N/C(/C(=O)O[C@H]1[C@@H](O[C@H]2O[C@@H]([C@](O)([C@@H](OC(=O)C(C)C)C)[C@@H](C2)OC)C)[C@@H](O)[C@H]([C@]3(O)C(O)=C(C(=N)C(C3)=O)C(=O)O)O[C@@H]1COC(=O)C)=C\\C '
               'NAME: Paulomycin B REASON: WRONGLY CLASSIFIED Contains one '
               'isothiocyanate group (R-N=C=S)\n'
               ' * SMILES: '
               'S=C=N[C@H]1[C@](C=C)(CC[C@H]2[C@H]1C=3C4=C(C=CC=C4)NC3C2(C)C)C '
               'NAME: 12-epi-fisherindole U isothiocyanate REASON: WRONGLY '
               'CLASSIFIED Contains one isothiocyanate group (R-N=C=S)\n'
               ' * SMILES: '
               'S=C=NCC1=CC=C(O[C@@H]2O[C@H]([C@H](OC(=O)C)[C@@H](O)[C@H]2O)C)C=C1 '
               'NAME: '
               "4-[(4'-o-acetyl-alpha-l-rhamnosyloxy)benzyl]isothiocyanate "
               'REASON: WRONGLY CLASSIFIED Contains one isothiocyanate group '
               '(R-N=C=S)\n'
               ' * SMILES: CCC(=O)N(C1CCN(CC1)CCC2=CC=C(C=C2)N=C=S)C3=CC=CC=C3 '
               'NAME: '
               'N-[1-[2-(4-isothiocyanatophenyl)ethyl]-4-piperidinyl]-N-phenylpropanamide '
               'REASON: WRONGLY CLASSIFIED Contains one isothiocyanate group '
               '(R-N=C=S)\n'
               ' * SMILES: '
               'S=C=N/C(/C(=O)O[C@H]1[C@@H](OC(=O)C)[C@@H](O)[C@H]([C@]2(O)C(O)=C(C(=N)C(C2)=O)C(=O)O)O[C@@H]1C)=C\\C '
               'NAME: Paulomycin G REASON: WRONGLY CLASSIFIED Contains one '
               'isothiocyanate group (R-N=C=S)\n'
               ' * SMILES: '
               'S=C=N/C(/C(=O)O[C@H]1[C@@H](O[C@@H]2O[C@H]([C@](O)([C@@H](OC(=O)C)C)[C@H](C2)OC)C)[C@@H](O)[C@H]([C@]3(O)C(O)=C(C(=N)C(C3)=O)C(=O)O)O[C@@H]1COC(=O)C)=C/C '
               'NAME: Paulomycin D REASON: WRONGLY CLASSIFIED Contains one '
               'isothiocyanate group (R-N=C=S)\n'
               'False negatives: NONE\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C1C=C(N=C2N1NN=N2)C3=CC=CC=C3',
                                     'name': '5-phenyl-1,7-dihydrotetrazolo[1,5-a]pyrimidine',
                                     'reason': 'No isothiocyanate group found'},
                                 {   'smiles': 'O(C(=O)CCCCCCCCCCC/C=C\\CCCCCCCC)[C@H](COC(=O)CCCCCCCCC/C=C\\CCCCCC)COC(=O)CCCCCCC/C=C\\CCCC',
                                     'name': 'TG(14:1(9Z)/22:1(13Z)/18:1(11Z))',
                                     'reason': 'No isothiocyanate group found'},
                                 {   'smiles': 'CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@@H](O[C@@H]3O[C@H](CO[C@H]4O[C@H](CO[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H](O)[C@H](O[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H]4O)[C@@H](O)[C@H](O[C@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]4O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4NC(C)=O)[C@@H]3O)[C@H](O)[C@H]2NC(C)=O)[C@@H]1O',
                                     'name': 'beta-D-GlcpNAc-(1->2)-alpha-D-Manp-(1->3)-{alpha-D-Manp-(1->3)-[alpha-D-Manp-(1->6)]-alpha-D-Manp-(1->6)}-beta-D-Manp-(1->4)-beta-GlcpNAc-(1->4)-beta-D-GlcpNAc',
                                     'reason': 'No isothiocyanate group found'},
                                 {   'smiles': '[C@@]12([C@@]([C@]([C@@H](CC1)C)(CCC(CC)C)C)(CCC[C@@H]2C)[H])C',
                                     'name': 'clerodane',
                                     'reason': 'No isothiocyanate group found'},
                                 {   'smiles': 'S(O[C@@H]1[C@H](O)[C@@H](O)[C@H](O[C@@H]([C@@H](O)[C@H](O)CO[C@]2(O[C@H]([C@H](NC(=O)C)[C@@H](O)C2)[C@H](O)[C@H](O)CO)C(O)=O)[C@@H](NC(=O)C)CO)O[C@@H]1CO)(O)(=O)=O',
                                     'name': '(2R,4S,5R,6R)-5-Acetamido-2-[(2R,3S,4R,5S)-5-acetamido-4-[(2R,3R,4R,5R,6R)-3,4-dihydroxy-6-(hydroxymethyl)-5-sulfooxyoxan-2-yl]oxy-2,3,6-trihydroxyhexoxy]-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'No isothiocyanate group found'},
                                 {   'smiles': 'O=C(NC(CC(C)C)C(=O)N[C@@H](CCCN=C(N)N)C=O)C(NC(=O)CC)CC(C)C',
                                     'name': 'Leupeptin Pr-LL',
                                     'reason': 'No isothiocyanate group found'},
                                 {   'smiles': 'COC1=CC=C(C=C1)CCNC(=O)C(=C2C3=CC=CC=C3C(=N2)NC(=O)C4=CC=CS4)C#N',
                                     'name': 'N-[3-[1-cyano-2-[2-(4-methoxyphenyl)ethylamino]-2-oxoethylidene]-1-isoindolyl]-2-thiophenecarboxamide',
                                     'reason': 'No isothiocyanate group found'},
                                 {   'smiles': 'O=C1O[C@@H](CC=2C1=C(OC)C(O)=CC2)CCC[C@@H](O)C',
                                     'name': 'Penicimarin C',
                                     'reason': 'No isothiocyanate group found'},
                                 {   'smiles': 'CO[C@@H]1[C@@H]2[C@@H](C[C@@H]3[C@@](O2)([C@H]([C@@H]([C@H](O3)C(=O)OC)O)O)O)OC1N4C=NC5=C4N=CN=C5N',
                                     'name': 'LSM-4497',
                                     'reason': 'No isothiocyanate group found'},
                                 {   'smiles': 'C[C@@H]1[C@H](O)CCC2=CC[C@H](C[C@]12C)C(C)=C',
                                     'name': '1-deoxycapsidiol',
                                     'reason': 'No isothiocyanate group '
                                               'found'}],
    'sample_false_negatives': [   {   'smiles': '[Cl-].CCN(CC)c1ccc2c(-c3ccc(cc3C(O)=O)N=C=S)c3ccc(cc3oc2c1)=[N+](CC)CC',
                                      'name': 'rhodamine B 5-isothiocyanate',
                                      'reason': 'No valid isothiocyanate '
                                                'groups found'},
                                  {   'smiles': 'CCN(CC)c1ccc2c(-c3cc(ccc3C(O)=O)N=C=S)c3ccc(cc3oc2c1)=[N+](CC)CC',
                                      'name': 'rhodamine B 6-isothiocyanate',
                                      'reason': 'No valid isothiocyanate '
                                                'groups found'},
                                  {   'smiles': 'Oc1c(Br)cc2c(Oc3c(Br)c(O)c(Br)cc3C22OC(=O)c3cc(ccc23)N=C=S)c1Br',
                                      'name': 'eosin 5-isothiocyanate',
                                      'reason': 'No valid isothiocyanate '
                                                'groups found'},
                                  {   'smiles': 'C1=CC(=CC=C1NC=2C=CC(N=C=S)=CC2)[N+](=O)[O-]',
                                      'name': 'amoscanate',
                                      'reason': 'Contains charged species that '
                                                'indicate different chemical '
                                                'class'}],
    'attempt': 2,
    'reasoning': 'Let me analyze the previous attempt and outcomes:\n'
                 '\n'
                 '1. The program correctly identified true positives but had '
                 'many false positives, indicating that our validation '
                 'criteria were too lenient.\n'
                 '\n'
                 '2. Looking at the false positives, many of them are complex '
                 'molecules with isothiocyanate groups but also other '
                 'structural features that might make them belong to different '
                 'chemical classes.\n'
                 '\n'
                 '3. The main issues appear to be:\n'
                 '   - Too simple validation of the R group - just checking if '
                 "it's connected to carbon isn't sufficient\n"
                 '   - No validation of the overall molecular properties\n'
                 '   - No checking for problematic features that would make '
                 'the molecule belong to a different class\n'
                 '\n'
                 '4. Improvements needed:\n'
                 '   - Add stricter validation of the isothiocyanate group '
                 'geometry\n'
                 '   - Better validation of the R group\n'
                 '   - Add checks for molecular properties and problematic '
                 'substructures\n'
                 '   - Ensure the isothiocyanate group is a primary feature of '
                 'the molecule\n'
                 '\n'
                 "Here's my improved version:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 27,
    'num_false_positives': 29,
    'num_true_negatives': 142240,
    'num_false_negatives': 4,
    'num_negatives': None,
    'precision': 0.48214285714285715,
    'recall': 0.8709677419354839,
    'f1': 0.6206896551724138,
    'accuracy': 0.9997680955727336,
    'negative_predictive_value': 0.9999718793059813}