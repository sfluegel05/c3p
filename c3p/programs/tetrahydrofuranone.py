"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
from rdkit import Chem

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone based on its SMILES string.
    A tetrahydrofuranone typically involves a five-membered lactone ring and
    may contain various levels of saturation and substitution.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is likely a tetrahydrofuranone derivative, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for various tetrahydrofuranone-like structures
    patterns = [Chem.MolFromSmarts(smarts) for smarts in [
        "O1C=CC(=O)C1",      # Unsaturated furanone with double bond
        "O1CC(=O)C=C1",      # Another variation with a double bond
        "O1CCOC1=O",         # Saturated with ester group
        "O1CCC(=O)C1",       # Simple saturated cyclic ester
        "O1C(=O)CC=C1",      # Alpha-beta unsaturated lactone
        "O1CC=C(O)C1",       # Unsaturated, with hydroxyl variations
        "O1C(=O)CCC1"        # Various saturated configurations
    ]]

    # Check if the molecule matches any of the defined tetrahydrofuranone patterns
    for pattern in patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains a tetrahydrofuranone-like structure"

    return False, "Does not contain a tetrahydrofuranone-like structure"

# Example usage
# Test the function with a sample SMILES string from known examples
test_smiles = "CCCCCCCC1OC(=O)C(=C)C1C(O)=O"
result, reason = is_tetrahydrofuranone(test_smiles)
print(f"Is Tetrahydrofuranone: {result}, Reason: {reason}")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:47016',
                          'name': 'tetrahydrofuranone',
                          'definition': 'Any oxolane having an oxo- '
                                        'substituent at any position on the '
                                        'tetrahydrofuran ring.',
                          'parents': ['CHEBI:26912'],
                          'xrefs': ['PMID:6047194'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.6222222222222222 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'C\\C=C(/C)C(=O)O[C@@H]1[C@@H]2OC[C@@]3(C)[C@H]2[C@](C)([C@@H](O)C[C@H]3OC(C)=O)[C@H]2CC[C@@]3(C)[C@H](C4COC(=O)C4)C(=O)C=C3[C@]12C '
               'NAME: 3-acetyl-7-tigloylvilasinin lactone REASON: CORRECT '
               'Contains a tetrahydrofuranone-like structure\n'
               ' * SMILES: '
               'CC(=O)CC[C@@H]1[C@H](OC(=O)C1=C)C1=C(C)[C@@H](O)CC1=O NAME: '
               'isosecotanapartholide REASON: CORRECT Contains a '
               'tetrahydrofuranone-like structure\n'
               ' * SMILES: CCCCC[C@H]1CCC(=O)O1 NAME: gamma-nonanolactone '
               'REASON: CORRECT Contains a tetrahydrofuranone-like structure\n'
               ' * SMILES: CC12CCCCC1C(=O)OC2=O NAME: methylhexahydrophthalic '
               'anhydride REASON: CORRECT Contains a tetrahydrofuranone-like '
               'structure\n'
               ' * SMILES: '
               '[H][C@@]1(C[C@@H](OC1=O)\\C=C\\CC\\C=C/C\\C=C/CCCCC)[C@](O)(CC(=O)O[C@H]([C@H](CC(O)=O)C(O)=O)C(O)=O)C(O)=O '
               'NAME: citrafungin B REASON: CORRECT Contains a '
               'tetrahydrofuranone-like structure\n'
               ' * SMILES: '
               '[C@]1([C@@](C(OC1)=O)([C@@H](CCCCC(C)C)O)[H])(CO)[H] NAME: '
               'SCB1 REASON: CORRECT Contains a tetrahydrofuranone-like '
               'structure\n'
               ' * SMILES: CCCCCCCCCCCC(=O)N[C@@H]1CCOC1=O NAME: '
               'N-dodecanoyl-L-homoserine lactone REASON: CORRECT Contains a '
               'tetrahydrofuranone-like structure\n'
               ' * SMILES: COc1ccc(cc1)S(=O)(=O)OCCC1CC(C)(C)OC1=O NAME: '
               '3-{[2-(4-methoxyphenylsulfonyl)oxy]ethyl}-5,5-dimethyldihydro-2(3H)-furanone '
               'REASON: CORRECT Contains a tetrahydrofuranone-like structure\n'
               ' * SMILES: O=C1OC(=O)C2CC=CCC12 NAME: '
               '1,2,3,6-tetrahydrophthalic anhydride REASON: CORRECT Contains '
               'a tetrahydrofuranone-like structure\n'
               ' * SMILES: N[C@H]1CCOC1=O NAME: L-homoserine lactone REASON: '
               'CORRECT Contains a tetrahydrofuranone-like structure\n'
               ' * SMILES: C(C=1C=C(C(=CC1)O)OS(O)(=O)=O)C2CCC(O2)=O NAME: '
               "5-(3',4'-dihydroxyphenyl)-gamma-valerolactone 3'-O-sulfate "
               'REASON: CORRECT Contains a tetrahydrofuranone-like structure\n'
               ' * SMILES: Oc1ccc(C=CC(=O)N[C@H]2CCOC2=O)cc1 NAME: '
               'N-(4-coumaroyl)-L-homoserine lactone REASON: CORRECT Contains '
               'a tetrahydrofuranone-like structure\n'
               ' * SMILES: '
               'C[C@@]12CC[C@H]3C(=CC[C@H]4C(C)(C)[C@H](O)CC[C@]34C)[C@@]1(C)CC[C@H]2C1COC(=O)C1 '
               'NAME: '
               '3-hydroxy-24,25,26,27-tetranortirucall-7-ene-23(21)-lactone '
               'REASON: CORRECT Contains a tetrahydrofuranone-like structure\n'
               ' * SMILES: C[C@@H]1CCC(=O)O1 NAME: (R)-gamma-valerolactone '
               'REASON: CORRECT Contains a tetrahydrofuranone-like structure\n'
               ' * SMILES: CC1(C)CC(CCl)C(=O)O1 NAME: '
               '3-chloromethyl-5,5-dimethylbutyrolactone REASON: CORRECT '
               'Contains a tetrahydrofuranone-like structure\n'
               ' * SMILES: OCC(O)C1OC(=O)C(=O)C1=O NAME: dehydroascorbic acid '
               'REASON: CORRECT Contains a tetrahydrofuranone-like structure\n'
               ' * SMILES: [C@]1([C@@](C(OC1)=O)([C@@H](CCCCCCC)O)[H])(CO)[H] '
               'NAME: SCB2 REASON: CORRECT Contains a tetrahydrofuranone-like '
               'structure\n'
               ' * SMILES: Cc1ccc(cc1)S(=O)(=O)OCC1CC(C)(C)OC1=O NAME: '
               '3-tosyloxymethyl-5,5-dimethylbutyrolactone REASON: CORRECT '
               'Contains a tetrahydrofuranone-like structure\n'
               ' * SMILES: O[C@@H]1COC(=O)[C@@H]1O NAME: erythrono-1,4-lactone '
               'REASON: CORRECT Contains a tetrahydrofuranone-like structure\n'
               ' * SMILES: CCCCCCCCC1OC(=O)C(=C)C1C(O)=O NAME: '
               '4-methylidene-2-octyl-5-oxotetrahydrofuran-3-carboxylic acid '
               'REASON: CORRECT Contains a tetrahydrofuranone-like structure\n'
               ' * SMILES: CC12CC=CCC1C(=O)OC2=O NAME: '
               'methyltetrahydrophthalic anhydride REASON: CORRECT Contains a '
               'tetrahydrofuranone-like structure\n'
               ' * SMILES: '
               'C([C@@H]1[C@H](COC1=O)CC2=CC=C3OCOC3=C2)C4=CC(=C(C(=C4)O)OC)OC '
               "NAME: (-)-5'-desmethylyatein REASON: CORRECT Contains a "
               'tetrahydrofuranone-like structure\n'
               ' * SMILES: COc1ccc(cc1)S(=O)(=O)OCC1CC(C)(C)OC1=O NAME: '
               '3-[(4-methoxyphenylsulfonyloxy)methyl]-5,5-dimethylbutyrolactone '
               'REASON: CORRECT Contains a tetrahydrofuranone-like structure\n'
               ' * SMILES: CC1(C)CC(CSC#N)C(=O)O1 NAME: '
               '3-thiocyanatomethyl-5,5-dimethylbutyrolactone REASON: CORRECT '
               'Contains a tetrahydrofuranone-like structure\n'
               ' * SMILES: O1C[C@@H](O)[C@H](O)C1=O NAME: D-threonolactone '
               'REASON: CORRECT Contains a tetrahydrofuranone-like structure\n'
               'False positives: SMILES: '
               'CC(C)[C@@]1(O)O[C@H]2C[C@@]3(C)OC(=CC3=O)\\C(C)=C/[C@H]3OC(=O)[C@@]1(C)[C@H]23 '
               'NAME: Eremantholide A REASON: WRONGLY CLASSIFIED Contains a '
               'tetrahydrofuranone-like structure\n'
               ' * SMILES: '
               'C[C@@H]1C[C@H](O)[C@@H]2[C@@H](CC[C@@H]3C[C@@H](O)CC(=O)O3)[C@@H](C)C=CC2=C1 '
               'NAME: monacolin J REASON: WRONGLY CLASSIFIED Contains a '
               'tetrahydrofuranone-like structure\n'
               ' * SMILES: '
               'COC(=O)CCc1c(C)c2=CC3=[N+]4C(=Cc5c(C)c(C=C)c6C=C7[N+]8=C(C=c1n2[Fe--]48n56)[C@]1(CCC(=O)O1)[C@@]7(C)O)C(C=C)=C3C '
               'NAME: cis-heme d hydroxychlorin gamma-spirolactone methyl '
               'ester REASON: WRONGLY CLASSIFIED Contains a '
               'tetrahydrofuranone-like structure\n'
               ' * SMILES: '
               'O=C1C2=C(C[C@]3(OC(=O)C=C3)OC2)C[C@H]4[C@]1(OC(=O)[C@@H]4C(O)CCCCCCC)C '
               'NAME: Sequoiamonascin A REASON: WRONGLY CLASSIFIED Contains a '
               'tetrahydrofuranone-like structure\n'
               ' * SMILES: '
               'O=C1[C@@]2(C(=O)O[C@]34[C@@H]2O[C@@H]([C@@H]1CO)[C@@H](C[C@@H](CCC[C@H]3C=C(C(=O)O)[C@H](C4)C)C)C)C '
               'NAME: Okilactomycin A REASON: WRONGLY CLASSIFIED Contains a '
               'tetrahydrofuranone-like structure\n'
               ' * SMILES: '
               'C[C@]12[C@H]3[C@H](C(O)=O)C4=C[C@H]5CC[C@@]4([C@@H](O)C5=C)[C@@]3(CC[C@@H]1O)OC2=O '
               'NAME: Antheridiogen-An REASON: WRONGLY CLASSIFIED Contains a '
               'tetrahydrofuranone-like structure\n'
               ' * SMILES: '
               '[H][C@@]1(CCC(=O)OC1)[C@@]1([H])CC[C@]2([H])[C@]3([H])CCC4CCCC[C@]4(C)[C@@]3([H])CC[C@]12C '
               'NAME: bufanolide REASON: WRONGLY CLASSIFIED Contains a '
               'tetrahydrofuranone-like structure\n'
               ' * SMILES: '
               'O=C1O[C@@H](C[C@H](C1)[C@@]2(C(=O)NC(C2=O)=C(C)C)C/C=C/CC(O)C)C '
               'NAME: Cladosporiumin C REASON: WRONGLY CLASSIFIED Contains a '
               'tetrahydrofuranone-like structure\n'
               ' * SMILES: '
               'ClC1=C(C(Cl)=C(O)C=C1O)/C(=C/C=C/[C@]2(OC(=O)C=C(C2)C)C)/C '
               'NAME: Cosmochlorin B REASON: WRONGLY CLASSIFIED Contains a '
               'tetrahydrofuranone-like structure\n'
               ' * SMILES: '
               'O=C1OC([C@@H]2CC[C@@]3([C@@H]([C@]2(C=C1)C)C[C@@]4(C(=O)[C@@]5([C@@]3(C(=O)O[C@H]5C)C4=C)O)C)C)(C)C '
               'NAME: Preaustinoid E REASON: WRONGLY CLASSIFIED Contains a '
               'tetrahydrofuranone-like structure\n'
               ' * SMILES: '
               'O=C1O[C@H]2C(C)(C)C([C@@]1([C@@]2(/C=C/C=C/C=C(\\[C@@H]3[C@H](C(=O)[C@]([C@H]3C)(O)C)C)/C)C)C)=O '
               'NAME: Emericelactone D REASON: WRONGLY CLASSIFIED Contains a '
               'tetrahydrofuranone-like structure\n'
               ' * SMILES: '
               'O([C@]1(C/C=C\\CCCCC)/C(/C(=O)C=C1)=C\\C=C/[C@@H]2OC(=O)CC2)C(=O)C '
               'NAME: clavulolactone I REASON: WRONGLY CLASSIFIED Contains a '
               'tetrahydrofuranone-like structure\n'
               ' * SMILES: '
               'O1[C@]23[C@@]4([C@@]5(C[C@@](OC(=O)C5=C)([C@@]3(OC(=O)[C@@]2(O)CC[C@]6([C@]([C@]1(O)C4=O)([C@H](O)C=C7[C@]6(C)C(=O)C=CC7)[H])[H])C)[H])C)[H] '
               'NAME: Physalin A REASON: WRONGLY CLASSIFIED Contains a '
               'tetrahydrofuranone-like structure\n'
               ' * SMILES: '
               'CO[C@H]1[C@H](O)[C@@H](CO)O[C@@H](O[C@H]2[C@H](O)CO[C@@H](O[C@@H]3[C@@H](C)O[C@@H](O[C@@H]4[C@@H](O)[C@@H](CO[C@H]4O[C@H]4CC[C@@]5(C)[C@@H](CC[C@@H]6C5=CC[C@]57[C@H](C(=O)C[C@@]65C)[C@](C)(CCCC(C)=C)OC7=O)C4(C)C)OS(O)(=O)=O)[C@H](O)[C@H]3O)[C@@H]2O)[C@@H]1O '
               'NAME: neothyonidioside REASON: WRONGLY CLASSIFIED Contains a '
               'tetrahydrofuranone-like structure\n'
               ' * SMILES: '
               'O=C1O[C@H]([C@H]([C@@H]2[C@@]3([C@@](C4=C([C@@]5([C@H](C([C@H](OC(=O)C)CC5)(C)C)CC4)C)CC3)(C)CC2)C)C)CC[C@H]1C '
               'NAME: Astraeusin Q REASON: WRONGLY CLASSIFIED Contains a '
               'tetrahydrofuranone-like structure\n'
               ' * SMILES: '
               'O=C1O[C@@H]([C@H]2C3(OC=4C=5C(C=CC4)=CC=CC5O3)CC[C@@H](C2)O)CC1 '
               'NAME: Anteaglonialide C REASON: WRONGLY CLASSIFIED Contains a '
               'tetrahydrofuranone-like structure\n'
               ' * SMILES: '
               'O1C(C2OC(CC2)C(O)CCCCCC(=O)CCCCC3OC(=O)C(C3)CC(=O)C)CCC1C(O)CCCCCCCCCC '
               'NAME: 9-Oxoasimicinone REASON: WRONGLY CLASSIFIED Contains a '
               'tetrahydrofuranone-like structure\n'
               ' * SMILES: '
               'O=C1OC(C2=CC[C@]3([C@H](C(=C2C1)C)C[C@](C(=O)O)(C)[C@@]45[C@@]3(C(=O)OC4)C(=O)[C@H](C)O5)C)(C)C '
               'NAME: '
               '(1R,2S,12S,14R,15R,17S)-2,6,6,11,14,17-hexamethyl-8,18,19-trioxo-7,16,20-trioxapentacyclo[13.3.3.01,15.02,12.05,10]henicosa-4,10-diene-14-carboxylic '
               'acid REASON: WRONGLY CLASSIFIED Contains a '
               'tetrahydrofuranone-like structure\n'
               ' * SMILES: '
               'O=C1O[C@@]2(O)[C@@H](O)C3=C[C@@](C=C)(C)CC=C3[C@@]45C2[C@]1(CC[C@@H]4C5)C '
               'NAME: Myrocin F REASON: WRONGLY CLASSIFIED Contains a '
               'tetrahydrofuranone-like structure\n'
               ' * SMILES: O=C1O[C@@H]([C@H](OC(=O)CC2=CC=CC=C2)C)CC1 NAME: '
               '1-(5-oxotetrahydrofuran-2-yl)ethyl 2-phenylacetate REASON: '
               'WRONGLY CLASSIFIED Contains a tetrahydrofuranone-like '
               'structure\n'
               ' * SMILES: '
               'CC(=O)O[C@H]1C[C@@](C)(O)[C@@H]2C[C@H]3[C@@H](OC(=O)C3=C)C=C(C)[C@@H]12 '
               'NAME: Gaillardin REASON: WRONGLY CLASSIFIED Contains a '
               'tetrahydrofuranone-like structure\n'
               ' * SMILES: '
               'O1C23C4(C(C5(C(CC4O)C(C(OC(=O)C)CC5OC(=O)C)(C)C)C)CCC2(C(OC(=O)C13)C=6C=COC6)C)C '
               'NAME: '
               '[13-Acetyloxy-7-(furan-3-yl)-19-hydroxy-1,8,12,16,16-pentamethyl-5-oxo-3,6-dioxapentacyclo[9.8.0.02,4.02,8.012,17]nonadecan-15-yl] '
               'acetate REASON: WRONGLY CLASSIFIED Contains a '
               'tetrahydrofuranone-like structure\n'
               ' * SMILES: '
               'COC(=O)[C@]12OC[C@@]34[C@H]1[C@@H](OC(=O)C=C(C)C)C(=O)O[C@@H]3C[C@H]1[C@H](C)C(=O)C(O[C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)=C[C@]1(C)[C@H]4[C@@H](O)[C@@H]2O '
               'NAME: Bruceoside A REASON: WRONGLY CLASSIFIED Contains a '
               'tetrahydrofuranone-like structure\n'
               ' * SMILES: O=C1O[C@@H](CC(=C1C)O)CCC(CC)C NAME: '
               'Alkyldihydropyrone C REASON: WRONGLY CLASSIFIED Contains a '
               'tetrahydrofuranone-like structure\n'
               ' * SMILES: '
               'O=C1OC(CC(C1)O)CCC2C3(C(OC(=O)C(C)C)CC(C)CC3C=CC2C)C NAME: '
               'Dihydromonacolin-MV REASON: WRONGLY CLASSIFIED Contains a '
               'tetrahydrofuranone-like structure\n'
               'False negatives: SMILES: C[C@@]1(O)OC[C@H](O)C1=O NAME: '
               '(2R,4S)-2-methyl-2,4-dihydroxydihydrofuran-3-one REASON: '
               'MISSED Does not contain a tetrahydrofuranone-like structure\n'
               ' * SMILES: '
               'C=1(C(=C(C(=C(C1Cl)C)C([H])=O)O)C/C=C(/CC/C=C(\\C)/[C@@]2(CC(C(O2)(C)C)=O)[H])\\C)O '
               'NAME: ascofuranone REASON: MISSED Does not contain a '
               'tetrahydrofuranone-like structure\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'P(OCC(OC(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCC/C=C\\C/C=C\\C/C=C\\CCCCCCCC)(OCCN(C)C)(O)=O',
                                     'name': 'Pe-nme2(20:3(5Z,8Z,11Z)/22:6(4Z,7Z,10Z,13Z,16Z,19Z))',
                                     'reason': 'Does not contain a '
                                               'tetrahydrofuranone-like '
                                               'structure'},
                                 {   'smiles': 'O1[C@]2([C@@H](O)[C@H](O)[C@@H](O)[C@@]1(OCC[C@H](CCC=C(C(OC[C@]3(O)[C@@H](O)[C@@](OC3)(OC2)[H])=O)C)C)[H])[H]',
                                     'name': 'Urceolide',
                                     'reason': 'Does not contain a '
                                               'tetrahydrofuranone-like '
                                               'structure'},
                                 {   'smiles': 'CC1=CC(=C(C=C1)C)C(=O)CSC2=NN=C(S2)C',
                                     'name': '1-(2,5-dimethylphenyl)-2-[(5-methyl-1,3,4-thiadiazol-2-yl)thio]ethanone',
                                     'reason': 'Does not contain a '
                                               'tetrahydrofuranone-like '
                                               'structure'},
                                 {   'smiles': 'O=C1C=C([C@@]2(C(C[C@@H](C2)O)(C)C)C)CC[C@@]1(O)C',
                                     'name': 'Enokipodin H',
                                     'reason': 'Does not contain a '
                                               'tetrahydrofuranone-like '
                                               'structure'},
                                 {   'smiles': 'CCCC[C@](Cn1cncn1)(C#N)c1ccc(Cl)cc1',
                                     'name': '(S)-myclobutanil',
                                     'reason': 'Does not contain a '
                                               'tetrahydrofuranone-like '
                                               'structure'},
                                 {   'smiles': 'O=C1O[C@@H]([C@H](NC(=O)[C@@H](NC(=O)CCCCC)CC(=O)O)C(=O)N[C@H](C(=O)N[C@H]2CC[C@H](N([C@H](C(N([C@H](C(N[C@H]1C(C)C)=O)CC3=CC=CC=C3)C)=O)CC(C)C)C2=O)O)CCCCN(C)C)C',
                                     'name': 'Cyanopeptolin D',
                                     'reason': 'Does not contain a '
                                               'tetrahydrofuranone-like '
                                               'structure'},
                                 {   'smiles': 'O1[C@]2(O)N([C@H](C(=O)N3[C@]2(CCC3)[H])C(C)C)C(=O)[C@@]1(NC(=O)[C@H]4CN([C@]5(C(=C4)C6=C7C(C5)=CNC7=CC=C6)[H])C)C',
                                     'name': 'Ergovaline',
                                     'reason': 'Does not contain a '
                                               'tetrahydrofuranone-like '
                                               'structure'},
                                 {   'smiles': 'OC(=O)CCCCCCC#CCCCC',
                                     'name': '8-tridecynoic acid',
                                     'reason': 'Does not contain a '
                                               'tetrahydrofuranone-like '
                                               'structure'},
                                 {   'smiles': 'CC(=O)N([O-])CCC[C@H](N)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CO)C(=O)N[C@H]([C@H](O)[C@H]1S[C@H]([C@H](O)[C@H]1O)n1ccc(=N)n(C)c1=O)C(O)=O',
                                     'name': 'desferrialbomycin epsilon(3-)',
                                     'reason': 'Does not contain a '
                                               'tetrahydrofuranone-like '
                                               'structure'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@H](O[C@@H]3[C@@H](CO)O[C@H](O[C@@H]4[C@@H](CO)O[C@H](O)[C@H](O)[C@H]4O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O',
                                     'name': 'alpha-maltotetraose',
                                     'reason': 'Does not contain a '
                                               'tetrahydrofuranone-like '
                                               'structure'}],
    'sample_false_negatives': [],
    'attempt': 3,
    'reasoning': 'To classify a chemical entity as a tetrahydrofuranone, we '
                 'need to ensure that the molecule features a five-membered '
                 'ring containing an oxygen atom and a carbonyl group '
                 '(lactone), which aligns with the cyclic ester structure '
                 'commonly known as a tetrahydrofuranone. In the previous '
                 'attempt, the pattern matching seems insufficiently '
                 'comprehensive in recognizing all variations of this '
                 'structure, which resulted in a high number of false '
                 'positives and some false negatives.\n'
                 '\n'
                 'In order to improve the SMILES pattern matching, we should:\n'
                 '\n'
                 '1. Consider both saturated and unsaturated furanone '
                 'structures.\n'
                 '2. Account for substitutions on the ring structure, which '
                 'could exist as side chains or ring variation.\n'
                 '3. Recognize the presence of furanone-like structures within '
                 'larger complex molecules.\n'
                 '\n'
                 'Here’s an updated approach with revised SMILES patterns that '
                 'aim to better capture the diversity of '
                 'tetrahydrofuranone-like structures:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 86,
    'num_false_positives': 100,
    'num_true_negatives': 5204,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.46236559139784944,
    'recall': 1.0,
    'f1': 0.6323529411764706,
    'accuracy': 0.9814471243042672,
    'negative_predictive_value': 1.0}