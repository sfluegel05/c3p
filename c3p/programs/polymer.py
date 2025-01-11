"""
Classifies: CHEBI:60027 polymer
"""
from rdkit import Chem

def is_polymer(smiles: str):
    """
    Determines if a molecule is a polymer based on its SMILES string.
    Polymers are characterized by repeating units (monomers), specific structural features,
    and often have long chain-like structures with possible branching or network formation.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is likely a polymer, False if unlikely.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a list of potential repeating unit SMARTS patterns
    repeating_unit_patterns = [
        Chem.MolFromSmarts("C(C(C)=O)"),  # ester unit -COOH (e.g., acrylates)
        Chem.MolFromSmarts("CNC(=O)"),   # amide unit -CONH- (e.g., nylons),
        Chem.MolFromSmarts("COC"),       # ether unit -O- (e.g., polyethylene oxide)
        Chem.MolFromSmarts("C(=C)C"),    # vinyl unit (e.g., polyvinyl chloride)
        Chem.MolFromSmarts("[*]([*])[*]") # generic branching pattern
    ]

    # Check for repeating units indicative of polymer structures
    for pattern in repeating_unit_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, f"Matches polymer-related structure: {Chem.MolToSmarts(pattern)}"

    # Calculate longest chain and check if it's indicative of a polymer
    chain_length_threshold = 30  # Adjusted threshold
    longest_chain = 0
    
    for atom in mol.GetAtoms():
        if atom.GetDegree() > 2:
            continue  # Skip obvious branch points for simple linear chains
        visited_atoms = set()
        to_visit = [(atom, 0)]  # Include branching depth

        while to_visit:
            current, depth = to_visit.pop()
            if current.GetIdx() in visited_atoms:
                continue
            visited_atoms.add(current.GetIdx())
            longest_chain = max(longest_chain, depth)
            for neighbor in current.GetNeighbors():
                if neighbor.GetIdx() not in visited_atoms:
                    if neighbor.GetDegree() > 2:
                        branch_depth = depth + 1
                    else:
                        branch_depth = depth + 1
                    to_visit.append((neighbor, branch_depth))

    if longest_chain >= chain_length_threshold:
        return True, f"Found long chain (length: {longest_chain}) typical for polymers"

    return False, "Does not fit polymer characteristics"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:60027',
                          'name': 'polymer',
                          'definition': 'A polymer is a mixture, which is '
                                        'composed of macromolecules of '
                                        'different kinds and which may be '
                                        'differentiated by composition, '
                                        'length, degree of branching etc..',
                          'parents': ['CHEBI:60004'],
                          'xrefs': ['Wikipedia:Polymer'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.2649006622516556 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'C=1(C(C[C@@H](CC1C)OC2O[C@@H]([C@H]([C@@H]([C@H]2O)O)O)COC(=O)CCCCCCCCCC(C)C)(C)C)/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C=3C(C[C@@H](CC3C)O)(C)C)\\C)\\C)/C)/C '
               'NAME: Thermozeaxanthin-13 REASON: CORRECT Matches '
               'polymer-related structure: COC\n'
               ' * SMILES: '
               'C(\\C=C(\\CC\\C=C(\\CCC=C(C)C)/C)/C)C\\C(\\C)=C/CC\\C(\\C)=C/CC\\C(\\C)=C/CC\\C(\\C)=C/CC\\C(\\C)=C/COP(OP(O)(O)=O)(O)=O '
               'NAME: ditrans,polycis-octaprenyl diphosphate REASON: CORRECT '
               'Found long chain (length: 37) typical for polymers\n'
               ' * SMILES: '
               'CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\Cc1ccccc1O '
               'NAME: 2-all-trans-decaprenylphenol REASON: CORRECT Found long '
               'chain (length: 45) typical for polymers\n'
               ' * SMILES: '
               'C(\\C=C(\\CC\\C=C(\\CCC=C(C)C)/C)/C)C\\C(\\C)=C/CC\\C(\\C)=C/CC\\C(\\C)=C/CC\\C(\\C)=C/COP(OP(O)(O)=O)(O)=O '
               'NAME: ditrans,polycis-heptaprenyl diphosphate REASON: CORRECT '
               'Found long chain (length: 33) typical for polymers\n'
               ' * SMILES: '
               '[C@@H](C(CO)(C)C)(C(NCCC([O-])=O)=O)O.[Ca+2].[C@@H](C(CO)(C)C)(C(NCCC([O-])=O)=O)O '
               'NAME: Calcium pantothenate REASON: CORRECT Matches '
               'polymer-related structure: CNC=O\n'
               ' * SMILES: '
               'C=1(C(C[C@@H](CC1C)OC2O[C@@H]([C@H]([C@@H]([C@H]2O)O)O)COC(=O)CCCCCCCCCCCC(C)C)(C)C)/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C=3C(CCCC3C)(C)C)\\C)\\C)/C)/C '
               'NAME: Thermocryptoxanthin-15 REASON: CORRECT Matches '
               'polymer-related structure: COC\n'
               ' * SMILES: '
               'N1(C(=O)C[C@@H](C([O-])=O)CC2=CC=CC=C2)C[C@]3([H])[C@@]([H])(CCCC3)C1.[Ca+2].N4(C(=O)C[C@@H](C([O-])=O)CC5=CC=CC=C5)C[C@]6([H])[C@@]([H])(CCCC6)C4 '
               'NAME: Kad 1229 REASON: CORRECT Matches polymer-related '
               'structure: CNC=O\n'
               ' * SMILES: '
               'CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\Cc1ccccc1O '
               'NAME: 2-octaprenylphenol REASON: CORRECT Found long chain '
               '(length: 37) typical for polymers\n'
               ' * SMILES: '
               'C1(O[C@@H]([C@H]([C@@H]([C@H]1O)O)O)COC(=O)CCCCCCCCCCCC(C)C)O[C@H]2CC(C(\\C=C\\C(=C\\C=C\\C(=C\\C=C\\C=C(\\C=C\\C=C(\\C=C\\C=3C(C[C@@H](CC3C)OC4O[C@@H]([C@H]([C@@H]([C@H]4O)O)O)COC(=O)CCCCCCCCCC(C)C)(C)C)/C)/C)\\C)\\C)=C(C2)C)(C)C '
               'NAME: Thermobiszeaxanthin-13-15 REASON: CORRECT Matches '
               'polymer-related structure: COC\n'
               ' * SMILES: '
               'C=1(C(C[C@@H](CC1C)OC2O[C@@H]([C@H]([C@@H]([C@H]2O)O)O)COC(=O)CCCCCCCCCC(C)C)(C)C)/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C=3C(CCCC3C)(C)C)\\C)\\C)/C)/C '
               'NAME: Thermocryptoxanthin-13 REASON: CORRECT Matches '
               'polymer-related structure: COC\n'
               ' * SMILES: '
               '[C@@H]1([C@@H]([C@H](OC([C@@H]1NC(C)=O)OP(OP(OC\\C=C(\\C)/CC\\C=C(\\C)/CC\\C=C(\\C)/CC\\C=C(\\C)/CC\\C=C(\\C)/CC\\C=C(\\C)/CC\\C=C(\\C)/CC\\C=C(\\C)/CC\\C=C(\\CC\\C=C(\\CCC=C(C)C)/C)/C)(O)=O)(O)=O)CO)O[C@H]2[C@@H]([C@H]([C@H](O)[C@H](O2)CO)O)NC(C)=O)OC(C(N[C@H](C(NC(CCC(N[C@H](C(N[C@@H](C(N[C@@H](C(O)=O)C)=O)C)=O)CCCCNC([C@@H](N)C)=O)=O)C(O)=O)=O)C)=O)C '
               'NAME: '
               'Undecaprenyl-diphospho-N-acetylmuramoyl-(N-acetylglucosamine)-L-alanyl-gamma-D-glutamyl-L-lysyl-(L-alanyl)-D-alanyl-D-alanine '
               'REASON: CORRECT Matches polymer-related structure: CNC=O\n'
               ' * SMILES: '
               'C(\\C=C(\\CC\\C=C(\\CCC=C(C)C)/C)/C)C\\C(\\C)=C/CC\\C(\\C)=C/CC\\C(\\C)=C/CC\\C(\\C)=C/CC\\C(\\C)=C/CC\\C(\\C)=C/COP(OP(O)(O)=O)(O)=O '
               'NAME: ditrans,polycis-nonaprenyl diphosphate REASON: CORRECT '
               'Found long chain (length: 41) typical for polymers\n'
               ' * SMILES: '
               '[C@@H]1([C@@H]([C@H](OC([C@@H]1NC(C)=O)OP(OP(OC\\C=C(\\C)/CC\\C=C(\\C)/CC\\C=C(\\C)/CC\\C=C(\\C)/CC\\C=C(\\C)/CC\\C=C(\\C)/CC\\C=C(\\C)/CC\\C=C(\\C)/CC\\C=C(\\CC\\C=C(\\CCC=C(C)C)/C)/C)(O)=O)(O)=O)CO)O[C@H]2[C@@H]([C@H]([C@H](O)[C@H](O2)CO)O)NC(C)=O)OC(C(N[C@H](C(NC(CCC(N[C@H](C(N[C@@H](C(N[C@@H](C(O)=O)C)=O)C)=O)CCCCNC([C@@H](NC([C@@H](N)C)=O)C)=O)=O)C(O)=O)=O)C)=O)C '
               'NAME: '
               'Undecaprenyl-diphospho-N-acetylmuramoyl-(N-acetylglucosamine)-L-alanyl-gamma-D-glutamyl-L-lysyl-(L-alanyl-L-alanyl)-D-alanyl-D-alanine '
               'REASON: CORRECT Matches polymer-related structure: CNC=O\n'
               ' * SMILES: '
               'C=1(C(C[C@@H](CC1C)OC2O[C@@H]([C@H]([C@@H]([C@H]2O)O)O)COC(=O)CCCCCCCCCCCC(C)C)(C)C)/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C=3C(C[C@@H](CC3C)O)(C)C)\\C)\\C)/C)/C '
               'NAME: Thermozeaxanthin-15 REASON: CORRECT Matches '
               'polymer-related structure: COC\n'
               ' * SMILES: '
               'CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\Cc1ccccc1O '
               'NAME: 2-all-trans-nonaprenylphenol REASON: CORRECT Found long '
               'chain (length: 41) typical for polymers\n'
               ' * SMILES: '
               '[C@@H]1(O[C@H]([C@H](O)[C@H]([C@H]1O)O)O[C@@H]([C@@H]([C@H](CO)O)O)[C@@H](CO)O)CO.O.O '
               'NAME: Lactitol dihydrate REASON: CORRECT Matches '
               'polymer-related structure: COC\n'
               ' * SMILES: '
               '[N+](O)(=O)[O-].C=1C=C2C(C(C=3[C@@](C[C@]4([H])[C@](C3O)(C(=O)C(=C(O)[C@H]4N(C)C)C(=O)NCN5CCCC5)O)([H])[C@]2(C)O)=O)=C(C1)O.O.[N+](O)(=O)[O-].C=6C=C7C(C(C=8[C@@](C[C@]9([H])[C@](C8O)(C(=O)C(=C(O)[C@H]9N(C)C)C(=O)NCN%10CCCC%10)O)([H])[C@]7(C)O)=O)=C(C6)O.O.O '
               'NAME: Rolitetracycline nitrate REASON: CORRECT Matches '
               'polymer-related structure: CC(C)=O\n'
               ' * SMILES: '
               'CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\Cc1ccccc1O '
               'NAME: 2-all-trans-heptaprenylphenol REASON: CORRECT Found long '
               'chain (length: 33) typical for polymers\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCC/C=C\\CCCCCCCCCCCCC/C=C/CCCCCCCCCCCCCCCCCC(C(C(=O)O)CCCCCCCCCCCCCCCCCCCCCC)O '
               'NAME: alpha-Semegma mycolic acid REASON: CORRECT Found long '
               'chain (length: 76) typical for polymers\n'
               ' * SMILES: '
               'C1(=C(C=C2C(=C1)C3=C(N2)[C@]4(C[C@@]5(C(=CO[C@H]([C@@]5(CN4CC3)[H])C)C(OCCN(C)C)=O)[H])[H])OC)OC.Cl.Cl '
               'NAME: Dimethylaminoethyl reserpilinate dihydrochloride REASON: '
               'CORRECT Matches polymer-related structure: COC\n'
               'False positives: SMILES: '
               'P(OCC(OC(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCC/C=C\\C/C=C\\C/C=C\\CCCCCCCC)(OCCN(C)C)(O)=O '
               'NAME: Pe-nme2(20:3(5Z,8Z,11Z)/22:6(4Z,7Z,10Z,13Z,16Z,19Z)) '
               'REASON: WRONGLY CLASSIFIED Matches polymer-related structure: '
               'COC\n'
               ' * SMILES: '
               'O1[C@]2([C@@H](O)[C@H](O)[C@@H](O)[C@@]1(OCC[C@H](CCC=C(C(OC[C@]3(O)[C@@H](O)[C@@](OC3)(OC2)[H])=O)C)C)[H])[H] '
               'NAME: Urceolide REASON: WRONGLY CLASSIFIED Matches '
               'polymer-related structure: COC\n'
               ' * SMILES: O=C1C=C([C@@]2(C(C[C@@H](C2)O)(C)C)C)CC[C@@]1(O)C '
               'NAME: Enokipodin H REASON: WRONGLY CLASSIFIED Matches '
               'polymer-related structure: CC(C)=O\n'
               ' * SMILES: '
               'O=C1O[C@@H]([C@H](NC(=O)[C@@H](NC(=O)CCCCC)CC(=O)O)C(=O)N[C@H](C(=O)N[C@H]2CC[C@H](N([C@H](C(N([C@H](C(N[C@H]1C(C)C)=O)CC3=CC=CC=C3)C)=O)CC(C)C)C2=O)O)CCCCN(C)C)C '
               'NAME: Cyanopeptolin D REASON: WRONGLY CLASSIFIED Matches '
               'polymer-related structure: CNC=O\n'
               ' * SMILES: '
               'O1[C@]2(O)N([C@H](C(=O)N3[C@]2(CCC3)[H])C(C)C)C(=O)[C@@]1(NC(=O)[C@H]4CN([C@]5(C(=C4)C6=C7C(C5)=CNC7=CC=C6)[H])C)C '
               'NAME: Ergovaline REASON: WRONGLY CLASSIFIED Matches '
               'polymer-related structure: CNC=O\n'
               ' * SMILES: '
               'CC(=O)N([O-])CCC[C@H](N)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CO)C(=O)N[C@H]([C@H](O)[C@H]1S[C@H]([C@H](O)[C@H]1O)n1ccc(=N)n(C)c1=O)C(O)=O '
               'NAME: desferrialbomycin epsilon(3-) REASON: WRONGLY CLASSIFIED '
               'Matches polymer-related structure: CNC=O\n'
               ' * SMILES: '
               'OC[C@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@H](O[C@@H]3[C@@H](CO)O[C@H](O[C@@H]4[C@@H](CO)O[C@H](O)[C@H](O)[C@H]4O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O '
               'NAME: alpha-maltotetraose REASON: WRONGLY CLASSIFIED Matches '
               'polymer-related structure: COC\n'
               ' * SMILES: '
               'SC[C@H](N)C(=O)N[C@H](C(=O)N[C@@H](CCC(=O)N)C(O)=O)CO NAME: '
               'Cys-Ser-Gln REASON: WRONGLY CLASSIFIED Matches polymer-related '
               'structure: CNC=O\n'
               ' * SMILES: '
               'C[C@@H]1CN([C@H](COC2=C(C=C(C=C2)NC(=O)NC3=CC=C(C=C3)Cl)C(=O)N(C[C@@H]1OC)C)C)C '
               'NAME: '
               '1-(4-chlorophenyl)-3-[(4S,7R,8R)-8-methoxy-4,5,7,10-tetramethyl-11-oxo-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]urea '
               'REASON: WRONGLY CLASSIFIED Matches polymer-related structure: '
               'CNC=O\n'
               ' * SMILES: '
               'C1C[C@@H]([C@@H](O[C@@H]1CCNC(=O)C2CCOCC2)CO)NC(=O)NC3=CC(=CC=C3)Cl '
               'NAME: '
               'N-[2-[(2S,5S,6R)-5-[[(3-chloroanilino)-oxomethyl]amino]-6-(hydroxymethyl)-2-oxanyl]ethyl]-4-oxanecarboxamide '
               'REASON: WRONGLY CLASSIFIED Matches polymer-related structure: '
               'CNC=O\n'
               ' * SMILES: CCCCN(CCCC)CCNC(=O)CNC(=O)C1=CC=CS1 NAME: '
               'N-[2-[2-(dibutylamino)ethylamino]-2-oxoethyl]-2-thiophenecarboxamide '
               'REASON: WRONGLY CLASSIFIED Matches polymer-related structure: '
               'CNC=O\n'
               ' * SMILES: '
               'C[C@@H](C1=CC=CC=C1)NC(=O)C[C@@H]2CC[C@H]([C@H](O2)CO)NC(=O)C3CCC3 '
               'NAME: '
               'N-[(2S,3R,6S)-2-(hydroxymethyl)-6-[2-oxo-2-[[(1S)-1-phenylethyl]amino]ethyl]-3-oxanyl]cyclobutanecarboxamide '
               'REASON: WRONGLY CLASSIFIED Matches polymer-related structure: '
               'CNC=O\n'
               ' * SMILES: '
               'OC1=CC=C2C3(OC(C4=C3C=C(C=C4)C(O)=O)=O)C=5C(OC2=C1)=CC(=CC5)O.OC1=CC=C2C3(OC(C4=C3C=CC(=C4)C(O)=O)=O)C=5C(OC2=C1)=CC(=CC5)O '
               'NAME: 5(6)-carboxyfluorescein REASON: WRONGLY CLASSIFIED '
               'Matches polymer-related structure: COC\n'
               ' * SMILES: '
               '[H][C@@]1(CC[C@@]2([H])[C@]3([H])CC(=O)[C@@]4([H])C[C@H](O)[C@H](O)C[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)[C@@H](O)[C@H](O)[C@@H](C)C(C)C '
               'NAME: castasterone REASON: WRONGLY CLASSIFIED Matches '
               'polymer-related structure: CC(C)=O\n'
               ' * SMILES: O=C1C=C([C@](O)(/C=C/C(/C(=O)O)=C/C)C(C1)(C)C)C '
               'NAME: (+)-(S)-1-methylabscisic-6-acid REASON: WRONGLY '
               'CLASSIFIED Matches polymer-related structure: CC(C)=O\n'
               ' * SMILES: '
               'CCOc1cc(ccc1Nc1ncc2N(C)C(=O)c3ccccc3N(C)c2n1)N1CCC(O)CC1 NAME: '
               'XMD8-92 REASON: WRONGLY CLASSIFIED Matches polymer-related '
               'structure: CNC=O\n'
               ' * SMILES: '
               'C[C@@H]1CN(C(=O)C2=C(N=CC(=C2)C#CC3CCCCC3)O[C@H]1CN(C)CC4=CN=CC=C4)[C@H](C)CO '
               'NAME: '
               '(2R,3R)-8-(2-cyclohexylethynyl)-5-[(2R)-1-hydroxypropan-2-yl]-3-methyl-2-[[methyl(3-pyridinylmethyl)amino]methyl]-3,4-dihydro-2H-pyrido[2,3-b][1,5]oxazocin-6-one '
               'REASON: WRONGLY CLASSIFIED Matches polymer-related structure: '
               'CNC=O\n'
               ' * SMILES: '
               'O1[C@@H](O[C@H]([C@@H](NC(=O)C)CO)[C@@H](O)[C@H](O)CO)[C@H](NC(=O)C)[C@@H](O)[C@H](O)[C@H]1CO '
               'NAME: '
               'N-[(2R,3R,4R,5S,6R)-2-[(2S,3R,4S,5R)-2-Acetamido-1,4,5,6-tetrahydroxyhexan-3-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]acetamide '
               'REASON: WRONGLY CLASSIFIED Matches polymer-related structure: '
               'CNC=O\n'
               ' * SMILES: '
               'S(C1=C(C=2C(=O)C3=C(C=C(C)C=C3O)OC2C=C1O)C(=O)OC)C4=C(C=5C(=O)C6=C(C=C(C)C=C6O)OC5C=C4O)C(=O)OC '
               'NAME: Castochrin REASON: WRONGLY CLASSIFIED Matches '
               'polymer-related structure: COC\n'
               ' * SMILES: '
               'S(OC1C2OC(=O)C=3C(C=4C(C(OCC1OC(O)C2O)=O)=CC(O)=C(O)C4O)=C(O)C(O)=C(O)C3)(O)(=O)=O '
               'NAME: '
               '{6,7,8,11,12,13,21,22-octahydroxy-3,16-dioxo-2,17,20-trioxatetracyclo[17.3.1.0?,?.0??,??]tricosa-4,6,8,10,12,14-hexaen-23-yl}oxidanesulfonic '
               'acid REASON: WRONGLY CLASSIFIED Matches polymer-related '
               'structure: COC\n'
               ' * SMILES: '
               'CC(C)[C@@]1(O)O[C@H]2C[C@@]3(C)OC(=CC3=O)\\C(C)=C/[C@H]3OC(=O)[C@@]1(C)[C@H]23 '
               'NAME: Eremantholide A REASON: WRONGLY CLASSIFIED Matches '
               'polymer-related structure: CC(C)=O\n'
               ' * SMILES: '
               'O1[C@@H]([C@@H](O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2NC(=O)C)CO)[C@H](O[C@H]3O[C@@H]([C@@H](O[C@@H]4O[C@@H]([C@@H](O)[C@H](O)[C@H]4NC(=O)C)CO)[C@H](O)[C@@H]3O)CO)[C@H](O)[C@@H]1O[C@H]5[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]5CO)O[C@H]6[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]6CO)O)CO[C@H]7O[C@@H]([C@@H](O)[C@H](O)[C@@H]7O[C@@H]8O[C@@H]([C@@H](O)[C@H](O)[C@H]8NC(=O)C)CO)CO '
               'NAME: '
               'N-[(2R,3R,4R,5S,6R)-5-[(2S,3R,4R,5S,6R)-3-Acetamido-5-[(2S,3S,4R,5R,6R)-5-[(2S,3R,4R,5S,6R)-3-acetamido-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-4-[(2R,3S,4R,5S,6R)-5-[(2S,3R,4R,5S,6R)-3-acetamido-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-3,4-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-6-[[(2S,3S,4S,5S,6R)-3-[(2S,3R,4R,5S,6R)-3-acetamido-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]-3-hydroxyoxan-2-yl]oxy-4-hydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-2,4-dihydroxy-6-(hydroxymethyl)oxan-3-yl]acetamide '
               'REASON: WRONGLY CLASSIFIED Matches polymer-related structure: '
               'CNC=O\n'
               ' * SMILES: '
               'C[C@H]1CN([C@@H](COC2=C(C=CC(=C2)NC(=O)C3CC3)C(=O)N(C[C@@H]1OC)C)C)C(=O)C4CCCCC4 '
               'NAME: '
               'N-[(5R,6S,9R)-8-[cyclohexyl(oxo)methyl]-5-methoxy-3,6,9-trimethyl-2-oxo-11-oxa-3,8-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]cyclopropanecarboxamide '
               'REASON: WRONGLY CLASSIFIED Matches polymer-related structure: '
               'CNC=O\n'
               ' * SMILES: '
               'C1[C@@H]2[C@H]([C@H](N2C(=O)CN1C(=O)C3=NC=CN=C3)CO)C4=CC=C(C=C4)C5=CC(=CC=C5)F '
               'NAME: '
               '(6S,7R,8S)-7-[4-(3-fluorophenyl)phenyl]-8-(hydroxymethyl)-4-[oxo(2-pyrazinyl)methyl]-1,4-diazabicyclo[4.2.0]octan-2-one '
               'REASON: WRONGLY CLASSIFIED Matches polymer-related structure: '
               'CNC=O\n'
               ' * SMILES: '
               '[NH2+]([C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(=O)N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(=O)N[C@H](C([O-])=C)CC(C)C)C(C)C)=O)CCCC[NH3+])=O)CC(C)C)=O)C(C)C)=O)CCC(N)=O)=O)CCCC[NH3+])=O)CCC(=O)N)=O)CCCC[NH3+])=O)C(C)O)=O)CCCC[NH3+])=O)CCCC[NH3+])=O)[C@H](CC)C)=O)C(C)C)=O)CCCNC(N)=[NH2+])=O)C(C)C)=O)CCCC[NH3+])=O)CC(C)C)=O)CO)=O)CCCC[NH3+])=O)C)=O)C)CCCC[NH3+])=O)C(C)O)=O)CC(C)C)=O)CC(=O)N)=O)CO)=O)C)C(=O)CCCC[C@@H]1[C@]2([C@@](CS1)(NC(N2)=O)[H])[H] '
               'NAME: HA Binder_2nd de novo REASON: WRONGLY CLASSIFIED Matches '
               'polymer-related structure: CNC=O\n'
               'False negatives: SMILES: '
               'CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\Cc1ccccc1O '
               'NAME: 2-hexaprenylphenol REASON: MISSED Does not fit polymer '
               'characteristics\n'
               ' * SMILES: '
               'C1CNCCN1.OC(CC(CC(=O)O)(O)C(=O)O)=O.C2CNCCN2.C3CNCCN3.OC(CC(CC(=O)O)(O)C(=O)O)=O '
               'NAME: Piperazine citrate REASON: MISSED Does not fit polymer '
               'characteristics\n'
               ' * SMILES: '
               'C1=CC(=C2C(=C1O)N=CC=C2)C(=O)C.S(=O)(=O)(O)O.C3=CC(=C4C(=C3O)N=CC=C4)C(=O)C '
               'NAME: Quinacetol sulfate REASON: MISSED Does not fit polymer '
               'characteristics\n'
               ' * SMILES: '
               '[As](=O)([O-])([O-])[O-].[Ca+2].[Ca+2].[Ca+2].[As](=O)([O-])([O-])[O-] '
               'NAME: Calcium arsenate REASON: MISSED Does not fit polymer '
               'characteristics\n'
               ' * SMILES: '
               'OC[C@H]([C@@H](O)[C@H](O)[C@@H](O)C([O-])=O)O.[Fe+3].OC[C@H]([C@@H](O)[C@H](O)[C@@H](O)C([O-])=O)O.OC[C@H]([C@@H](O)[C@H](O)[C@@H](O)C([O-])=O)O '
               'NAME: Ferric gluconate REASON: MISSED Does not fit polymer '
               'characteristics\n'
               ' * SMILES: '
               'C=1(N(C(N=C(N1)N)(C)C)C=2C=CC(=CC2)Cl)N.C3=CC=CC4=C3C=C(C(=C4CC=5C6=C(C=C(C5O)C(=O)O)C=CC=C6)O)C(=O)O.C=7(N(C(N=C(N7)N)(C)C)C=8C=CC(=CC8)Cl)N '
               'NAME: Cycloguanil pamoate REASON: MISSED Does not fit polymer '
               'characteristics\n'
               ' * SMILES: '
               'C(\\C=C(\\CC\\C=C(\\CCC=C(C)C)/C)/C)C\\C(\\C)=C/CC\\C(\\C)=C/CC\\C(\\C)=C/COP(OP(O)(O)=O)(O)=O '
               'NAME: ditrans,polycis-hexaprenyl diphosphate REASON: MISSED '
               'Does not fit polymer characteristics\n'
               ' * SMILES: '
               'C1CC=2C(N(C=3C1=CC(=CC3)Cl)CCCN4CCC(CC4)(N5CCCCC5)C(N)=O)=CC=CC2.Cl.Cl '
               'NAME: Clocapramine dihydrochloride REASON: MISSED Does not fit '
               'polymer characteristics\n'
               ' * SMILES: '
               'C1=C(C=CC(=C1[N+](=O)[O-])C(NC2=CC=C(C=C2)C=3NCCN3)=O)C(NC4=CC=C(C=C4)C=5NCCN5)=O.Cl.Cl '
               'NAME: HR1917 REASON: MISSED Does not fit polymer '
               'characteristics\n'
               ' * SMILES: '
               'C=1C=C2C(N(C3=C2CN(C)CC3)CC4=CC=CC=C4)=CC1.C5=CC=C6C(=C5S(O)(=O)=O)C=CC=C6S(O)(=O)=O.C=7C=C8C(N(C9=C8CN(C)CC9)CC%10=CC=CC=C%10)=CC7 '
               'NAME: Mebhydrolin napadisilate REASON: MISSED Does not fit '
               'polymer characteristics\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [],
    'sample_false_negatives': [],
    'attempt': 3,
    'reasoning': 'To improve the classification of polymers based on previous '
                 'attempts and outcomes, we should focus on refining our '
                 "method for identifying polymer characteristics. Here's a "
                 'breakdown of our observations and potential improvements:\n'
                 '\n'
                 '1. **Revisiting Polymer Definitions**: Polymers are '
                 'typically characterized by repeating structures (monomers) '
                 'and often form long chains. Our current method checks for '
                 'certain functional groups and long chains, which is a good '
                 'starting point but might not be exhaustive for all polymer '
                 'types.\n'
                 '\n'
                 '2. **Improving Repeating Unit Detection**: We used a '
                 'predefined set of repeating units commonly found in '
                 'polymers. However, we should expand this list to include a '
                 'wider variety of monomer patterns, possibly incorporating '
                 'additional patterns for common polymers like polyesters, '
                 'polyamides, and polyethers.\n'
                 '\n'
                 '3. **Chain Length Adjustment**: The chain length threshold '
                 'of 30 might not always capture all polymers, especially if '
                 'they are branched or have more complex configurations. '
                 'Adjust this threshold dynamically based on the presence of '
                 'branching or cross-linking elements.\n'
                 '\n'
                 '4. **Branch Structure Recognition**: Polymers can have '
                 'complex structures including branches or network shapes. '
                 'Incorporate additional checks for branch points and '
                 'network-connect structures, which might be characteristic of '
                 'certain polymer classes.\n'
                 '\n'
                 "5. **Structural Flexibility**: Utilize RDKit's capabilities "
                 'to leverage molecular descriptors or fingerprints to better '
                 'capture and quantify the complexity associated with polymer '
                 'structures.\n'
                 '\n'
                 "Based on these insights, here's an enhanced version of the "
                 '`is_polymer` function that attempts to address these points:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 31,
    'num_false_positives': 100,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.2366412213740458,
    'recall': 1.0,
    'f1': 0.38271604938271603,
    'accuracy': 0.2366412213740458,
    'negative_predictive_value': 0.0}