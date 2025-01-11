"""
Classifies: CHEBI:26167 polar amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

def is_polar_amino_acid(smiles: str):
    """
    Determines if a molecule is a polar amino acid based on its SMILES string.
    Polar amino acids have side chains capable of forming hydrogen bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polar amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check molecular weight - amino acids should be relatively small
    mol_weight = Descriptors.ExactMolWt(mol)
    if mol_weight > 250:  # Most amino acids are under 200 Da
        return False, "Molecule too large to be a single amino acid"

    # Check for basic amino acid structure (NH2-CH-COOH)
    amino_acid_pattern = Chem.MolFromSmarts("[NX3,NX4+][CX4][CX3](=[OX1])[OX2H,OX1-]")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "No amino acid backbone found"

    # Check number of amino groups and carboxyl groups
    amino_pattern = Chem.MolFromSmarts("[NX3,NX4+][CX4][CX3](=[OX1])[OX2H,OX1-]")
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H,OX1-]")
    
    amino_count = len(mol.GetSubstructMatches(amino_pattern))
    carboxyl_count = len(mol.GetSubstructMatches(carboxyl_pattern))
    
    if amino_count > 1 or carboxyl_count > 2:  # Allow up to 2 carboxyls for asp/glu
        return False, "Multiple amino acid residues detected - likely a peptide"

    # Initialize patterns for polar groups in side chains
    patterns = {
        "hydroxyl": "[OX2H]",  # -OH (serine, threonine, tyrosine)
        "amide": "[NX3H2][CX3](=[OX1])",  # -CONH2 (asparagine, glutamine)
        "carboxyl": "[CX3](=[OX1])[OX2H,OX1-]",  # -COOH (aspartic acid, glutamic acid)
        "basic_N": "[NX3;H2,H1;!$(NC=O)]",  # Basic N (lysine, arginine)
        "guanidino": "[NX3][CX3](=[NX2])[NX3]",  # Guanidino group (arginine)
        "imidazole": "c1c[nH]cn1",  # Imidazole ring (histidine)
        "thiol": "[SX2H]",  # -SH (cysteine)
        "indole": "c1ccc2c(c1)[nH]cc2",  # Indole ring (tryptophan)
        "phenol": "c1cc(O)ccc1"  # Phenol group (tyrosine)
    }

    # Get backbone atoms to exclude them from side chain search
    backbone_matches = mol.GetSubstructMatches(amino_acid_pattern)
    if not backbone_matches:
        return False, "No amino acid backbone found"
    
    backbone_atoms = set(backbone_matches[0])
    
    # Check each polar pattern
    found_polar_groups = []
    for group_name, pattern in patterns.items():
        pattern_mol = Chem.MolFromSmarts(pattern)
        matches = mol.GetSubstructMatches(pattern_mol)
        
        # Only count matches that aren't part of the backbone
        side_chain_matches = [match for match in matches if not all(atom_idx in backbone_atoms for atom_idx in match)]
        
        if side_chain_matches:
            found_polar_groups.append(group_name)

    if not found_polar_groups:
        # Special case for tryptophan - check for indole NH
        if mol.HasSubstructMatch(Chem.MolFromSmarts("c1ccc2c(c1)[nH]cc2")):
            return True, "Polar amino acid with indole NH group capable of hydrogen bonding"
        return False, "No polar side chain groups found"

    # Construct reason string
    polar_groups_str = ", ".join(found_polar_groups)
    return True, f"Polar amino acid with {polar_groups_str} group(s) in side chain"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26167',
                          'name': 'polar amino acid',
                          'definition': 'Any amino acid whose side chain is '
                                        'capable of forming one or more '
                                        'hydrogen bonds.',
                          'parents': ['CHEBI:33709'],
                          'xrefs': [   'MetaCyc:Polar-amino-acids',
                                       'PMID:12016058'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.4069767441860465 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: CC(CCNC(N)=N)C(N)C(O)=O NAME: '
               '3-methylarginine REASON: CORRECT Polar amino acid with '
               'basic_N, guanidino group(s) in side chain\n'
               ' * SMILES: NCCCC[C@@H](N)C(O)=O NAME: D-lysine REASON: CORRECT '
               'Polar amino acid with basic_N group(s) in side chain\n'
               ' * SMILES: '
               'OC(=O)[C@@H](N)CC1=C(C(=C(O)C(=C1[2H])[2H])[2H])[2H] NAME: '
               'L-tyrosine-d4 REASON: CORRECT Polar amino acid with hydroxyl '
               'group(s) in side chain\n'
               ' * SMILES: NC(CCC(O)=O)C(O)=O NAME: glutamic acid REASON: '
               'CORRECT Polar amino acid with hydroxyl, carboxyl group(s) in '
               'side chain\n'
               ' * SMILES: N[C@@H](CC(N)=O)C(O)=O NAME: L-asparagine REASON: '
               'CORRECT Polar amino acid with amide group(s) in side chain\n'
               ' * SMILES: N[C@@H](CCCNC(N)=N)C(O)=O NAME: L-arginine REASON: '
               'CORRECT Polar amino acid with basic_N, guanidino group(s) in '
               'side chain\n'
               ' * SMILES: NC(Cc1c[nH]cn1)C(O)=O NAME: histidine REASON: '
               'CORRECT Polar amino acid with imidazole group(s) in side '
               'chain\n'
               ' * SMILES: '
               'O(C(=O)[C@@](N([2H])[2H])(C(C(O[2H])=O)([2H])[2H])[2H])[2H] '
               'NAME: L-aspartic acid-d7 REASON: CORRECT Polar amino acid with '
               'hydroxyl, carboxyl group(s) in side chain\n'
               ' * SMILES: N[C@H](CO)C(O)=O NAME: D-serine REASON: CORRECT '
               'Polar amino acid with hydroxyl group(s) in side chain\n'
               ' * SMILES: N[C@H](Cc1c[nH]cn1)C(O)=O NAME: D-histidine REASON: '
               'CORRECT Polar amino acid with imidazole group(s) in side '
               'chain\n'
               ' * SMILES: N[C@@H](CC(O)=O)C(O)=O NAME: L-aspartic acid '
               'REASON: CORRECT Polar amino acid with hydroxyl, carboxyl '
               'group(s) in side chain\n'
               ' * SMILES: N[C@@H](Cc1c[nH]cn1)C(O)=O NAME: L-histidine '
               'REASON: CORRECT Polar amino acid with imidazole group(s) in '
               'side chain\n'
               ' * SMILES: '
               'OC(=O)[C@](N([2H])[2H])(C(C(CN=C(N)N)([2H])[2H])([2H])[2H])[2H] '
               'NAME: L-arginine-d7 REASON: CORRECT Polar amino acid with '
               'basic_N, guanidino group(s) in side chain\n'
               ' * SMILES: N[C@H](CC(O)=O)C(O)=O NAME: D-aspartic acid REASON: '
               'CORRECT Polar amino acid with hydroxyl, carboxyl group(s) in '
               'side chain\n'
               ' * SMILES: NC(C[13C](O)=O)[13C](O)=O NAME: aspartic '
               'acid-1,4-(13)C2 REASON: CORRECT Polar amino acid with '
               'hydroxyl, carboxyl group(s) in side chain\n'
               ' * SMILES: N[C@H](CS)C(O)=O NAME: D-cysteine REASON: CORRECT '
               'Polar amino acid with thiol group(s) in side chain\n'
               ' * SMILES: C[C@@H](O)[C@H](N)C(O)=O NAME: L-threonine REASON: '
               'CORRECT Polar amino acid with hydroxyl group(s) in side chain\n'
               ' * SMILES: NC(CO)C(O)=O NAME: serine REASON: CORRECT Polar '
               'amino acid with hydroxyl group(s) in side chain\n'
               ' * SMILES: N[C@@H](CCC(N)=O)C(O)=O NAME: L-glutamine REASON: '
               'CORRECT Polar amino acid with amide group(s) in side chain\n'
               ' * SMILES: NC(CC(O)=O)C(O)=O NAME: aspartic acid REASON: '
               'CORRECT Polar amino acid with hydroxyl, carboxyl group(s) in '
               'side chain\n'
               ' * SMILES: N[C@H](C[13C](O)=O)[13C](O)=O NAME: D-aspartic '
               'acid-1,4-(13)C2 REASON: CORRECT Polar amino acid with '
               'hydroxyl, carboxyl group(s) in side chain\n'
               ' * SMILES: N[C@@H](CCC(O)=O)C(O)=O NAME: L-glutamic acid '
               'REASON: CORRECT Polar amino acid with hydroxyl, carboxyl '
               'group(s) in side chain\n'
               ' * SMILES: NC(Cc1ccc(O)cc1)C(O)=O NAME: tyrosine REASON: '
               'CORRECT Polar amino acid with hydroxyl group(s) in side chain\n'
               ' * SMILES: C[C@H](O)[C@@H](N)C(O)=O NAME: D-threonine REASON: '
               'CORRECT Polar amino acid with hydroxyl group(s) in side chain\n'
               ' * SMILES: N[C@@H](C[13C](O)=O)[13C](O)=O NAME: L-aspartic '
               'acid-1,4-(13)C2 REASON: CORRECT Polar amino acid with '
               'hydroxyl, carboxyl group(s) in side chain\n'
               'False positives: SMILES: '
               'O[C@@H]([C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)N[C@@H](CC(C)C)C(O)=O)C '
               'NAME: Glu-Thr-Leu REASON: WRONGLY CLASSIFIED Polar amino acid '
               'with hydroxyl, carboxyl, basic_N group(s) in side chain\n'
               ' * SMILES: '
               'O=C(N[C@@H](CCC(O)=O)C(O)=O)[C@@H](NC(=O)[C@@H](N)CCCN=C(N)N)CC(C)C '
               'NAME: Arg-Leu-Glu REASON: WRONGLY CLASSIFIED Polar amino acid '
               'with hydroxyl, carboxyl, basic_N, guanidino group(s) in side '
               'chain\n'
               ' * SMILES: '
               'O=C(N[C@@H](C(C)C)C(O)=O)[C@@H](NC(=O)[C@@H](N)CC1=CC=C(O)C=C1)C '
               'NAME: Tyr-Ala-Val REASON: WRONGLY CLASSIFIED Polar amino acid '
               'with hydroxyl, basic_N group(s) in side chain\n'
               ' * SMILES: '
               'S(CC[C@H](N)C(=O)N[C@@H]([C@H](O)C)C(=O)N[C@@H](CC1=CC=CC=C1)C(O)=O)C '
               'NAME: Met-Thr-Phe REASON: WRONGLY CLASSIFIED Polar amino acid '
               'with hydroxyl, basic_N group(s) in side chain\n'
               ' * SMILES: '
               'O=C(N[C@@H](CC(O)=O)C(O)=O)[C@@H](NC(=O)[C@@H](N)CC1=CC=C(O)C=C1)CC=2NC=NC2 '
               'NAME: Tyr-His-Asp REASON: WRONGLY CLASSIFIED Polar amino acid '
               'with hydroxyl, carboxyl, basic_N, imidazole group(s) in side '
               'chain\n'
               ' * SMILES: '
               'O[C@@H]([C@H](NC(=O)[C@@H](N)CC=1C=2C(NC1)=CC=CC2)C(=O)N[C@@H](C)C(O)=O)C '
               'NAME: Trp-Thr-Ala REASON: WRONGLY CLASSIFIED Polar amino acid '
               'with hydroxyl, basic_N group(s) in side chain\n'
               ' * SMILES: '
               'SC[C@H](N)C(=O)NCC(=O)N[C@@H](CC=1C=2C(NC1)=CC=CC2)C(O)=O '
               'NAME: Cys-Gly-Trp REASON: WRONGLY CLASSIFIED Polar amino acid '
               'with basic_N, thiol group(s) in side chain\n'
               ' * SMILES: '
               'SC[C@H](NC(=O)[C@@H](N)[C@H](O)C)C(=O)N1[C@H](CCC1)C(O)=O '
               'NAME: Thr-Cys-Pro REASON: WRONGLY CLASSIFIED Polar amino acid '
               'with hydroxyl, basic_N, thiol group(s) in side chain\n'
               ' * SMILES: '
               'SC[C@H](N)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H]([C@H](CC)C)C(O)=O '
               'NAME: Cys-Leu-Ile REASON: WRONGLY CLASSIFIED Polar amino acid '
               'with basic_N, thiol group(s) in side chain\n'
               ' * SMILES: CC(C)C[C@H](N)C(=O)NCC(O)=O NAME: Leu-Gly REASON: '
               'WRONGLY CLASSIFIED Polar amino acid with basic_N group(s) in '
               'side chain\n'
               ' * SMILES: '
               'ClC(Cl)CCCCC(N)C(O)C(=O)N[C@H](C(=O)N([C@H](C(=O)N1[C@H](C(=O)O)CCC1)CC2=CC=C(O)C=C2)C)CC3=CC=C(O)C=C3 '
               'NAME: Microginin 680 REASON: WRONGLY CLASSIFIED Polar amino '
               'acid with hydroxyl, basic_N group(s) in side chain\n'
               ' * SMILES: CC(=O)NC(CCC(N)=O)C(O)=O NAME: N(2)-acetylglutamine '
               'REASON: WRONGLY CLASSIFIED Polar amino acid with amide '
               'group(s) in side chain\n'
               ' * SMILES: '
               'O=C(N[C@@H](CC1=CC=C(O)C=C1)C(O)=O)[C@@H](NC(=O)[C@@H](N)CC(=O)N)CCC(=O)N '
               'NAME: Asn-Gln-Tyr REASON: WRONGLY CLASSIFIED Polar amino acid '
               'with hydroxyl, amide, basic_N group(s) in side chain\n'
               ' * SMILES: '
               'O=C(N[C@@H](CC(O)=O)C(O)=O)[C@@H](NC(=O)[C@@H](N)CCC(=O)N)CC=1C=2C(NC1)=CC=CC2 '
               'NAME: Gln-Trp-Asp REASON: WRONGLY CLASSIFIED Polar amino acid '
               'with hydroxyl, amide, carboxyl, basic_N group(s) in side '
               'chain\n'
               ' * SMILES: '
               'S(CC[C@H](NC(=O)[C@@H](NC(=O)[C@@H](N)CO)C(C)C)C(O)=O)C NAME: '
               'Ser-Val-Met REASON: WRONGLY CLASSIFIED Polar amino acid with '
               'hydroxyl, basic_N group(s) in side chain\n'
               ' * SMILES: '
               'P(OC[C@H](OC(=O)CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCCCCCCCC/C=C\\CCCCCCCC)(OC[C@H](N)C(O)=O)(O)=O '
               'NAME: PS(20:1(11Z)/18:3(9Z,12Z,15Z)) REASON: WRONGLY '
               'CLASSIFIED Polar amino acid with hydroxyl group(s) in side '
               'chain\n'
               ' * SMILES: '
               'O=C(N1[C@@H](CCC1)C(=O)N[C@@H](C)C(O)=O)[C@@H](N)CC2=CC=C(O)C=C2 '
               'NAME: Tyr-Pro-Ala REASON: WRONGLY CLASSIFIED Polar amino acid '
               'with hydroxyl, basic_N group(s) in side chain\n'
               ' * SMILES: '
               '[C@H](NC([C@H](CC([O-])=O)[NH3+])=O)(CCCNC(N)=[NH2+])C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N1[C@@H](CCC1)C(N[C@@H](CC2=CC=CC=C2)C(N[C@@H](CC3=CNC=N3)C(N[C@H](C([O-])=O)CC(C)C)=O)=O)=O)=O)CC4=CN=CN4)=O)[C@H](CC)C)=O)CC5=CC=C(C=C5)O)=O)C(C)C)=O '
               'NAME: angiotensin I dizwitterion REASON: WRONGLY CLASSIFIED '
               'Polar amino acid with hydroxyl, carboxyl, basic_N, imidazole '
               'group(s) in side chain\n'
               ' * SMILES: '
               'CCCCCCCC\\C=C/CCCCCCCC(=O)OC[C@H](COP(O)(=O)OC[C@H](N)C(O)=O)OC(=O)CCCCCCC\\C=C/C\\C=C/CCCCC '
               'NAME: 1-oleoyl-2-linoleoyl-sn-glycero-3-phospho-L-serine '
               'REASON: WRONGLY CLASSIFIED Polar amino acid with hydroxyl '
               'group(s) in side chain\n'
               ' * SMILES: '
               'O=C(N[C@@H]([C@H](CC)C)C(O)=O)[C@@H](NC(=O)[C@@H](N)CCC(O)=O)CC(=O)N '
               'NAME: Glu-Asn-Ile REASON: WRONGLY CLASSIFIED Polar amino acid '
               'with hydroxyl, amide, carboxyl, basic_N group(s) in side '
               'chain\n'
               ' * SMILES: '
               '[NH3+][C@@H](CCSC[C@H]1OC(O)[C@H](O)[C@@H]1O)C([O-])=O NAME: '
               'S-(5-deoxy-D-ribos-5-yl)-L-homocysteine zwitterion REASON: '
               'WRONGLY CLASSIFIED Polar amino acid with hydroxyl group(s) in '
               'side chain\n'
               ' * SMILES: '
               'O=C(N[C@@H](CCC(=O)N)C(O)=O)[C@H]1N(CCC1)C(=O)[C@@H](N)CO '
               'NAME: Ser-Pro-Gln REASON: WRONGLY CLASSIFIED Polar amino acid '
               'with hydroxyl, amide, basic_N group(s) in side chain\n'
               ' * SMILES: '
               'SC[C@H](NC(=O)[C@@H](N)CC(O)=O)C(=O)N[C@@H](CC(=O)N)C(O)=O '
               'NAME: Asp-Cys-Asn REASON: WRONGLY CLASSIFIED Polar amino acid '
               'with hydroxyl, amide, carboxyl, basic_N, thiol group(s) in '
               'side chain\n'
               ' * SMILES: '
               '[C@@H]1([C@@H]([C@H]([C@@H]([C@H](O1)CO)O)O)NC(=O)C)O[C@H]2[C@@H]([C@H](C(O[C@@H]2CO)O)NC(=O)C)O[C@@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)O)CCC(N[C@H](C(O)=O)CCCC(N)C(=O)O)=O)C)C '
               'NAME: '
               'N-acetyl-D-glucosaminyl-(1->4)-N-acetylmuramoyl-L-alanyl-D-glutamyl-6-carboxy-L-lysine '
               'REASON: WRONGLY CLASSIFIED Polar amino acid with hydroxyl '
               'group(s) in side chain\n'
               ' * SMILES: '
               'O=C(N[C@@H](CCCCN)C(O)=O)[C@@H](NC(=O)[C@@H](N)[C@H](O)C)CC1=CC=CC=C1 '
               'NAME: Thr-Phe-Lys REASON: WRONGLY CLASSIFIED Polar amino acid '
               'with hydroxyl, basic_N group(s) in side chain\n'
               'False negatives: SMILES: N[C@H](Cc1c[nH]c2ccccc12)C(O)=O NAME: '
               'D-tryptophan REASON: MISSED No polar side chain groups found\n'
               ' * SMILES: NC(Cc1c[nH]c2ccccc12)C(O)=O NAME: tryptophan '
               'REASON: MISSED No polar side chain groups found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'Cc1cc(C)c2Cc3c(C)cc(C)c(Cc4c(C)cc(C)c(Cc5c(C)cc(C)c(Cc1c2C)c5C)c4C)c3C',
                                     'name': '4,6,10,12,16,18,22,24,25,26,27,28-dodecamethylcalix[4]arene',
                                     'reason': 'Molecule too large to be a '
                                               'single amino acid'},
                                 {   'smiles': 'O=C1C(=C2C=C3[C@]([C@@H](C(C)C)[C@@H]([C@H]3O)OC(=O)C)(C)CC[C@]2(C)CC1)COC(=O)C',
                                     'name': 'Dahliane E',
                                     'reason': 'Molecule too large to be a '
                                               'single amino acid'},
                                 {   'smiles': 'O[C@@H]([C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)N[C@@H](CC(C)C)C(O)=O)C',
                                     'name': 'Glu-Thr-Leu',
                                     'reason': 'Molecule too large to be a '
                                               'single amino acid'},
                                 {   'smiles': 'CCOc1ccc(NC(=O)C(C)O)cc1',
                                     'name': 'p-Lactophenetide',
                                     'reason': 'No amino acid backbone found'},
                                 {   'smiles': 'O=C1NCC=2C1=C3C(N([C@@H]4O[C@H]([C@@H](O)[C@H]([C@H]4OC)O)C)C5=C3C=CC=C5)=C6NC7=C(C26)C=CC=C7',
                                     'name': "3'-epi-5'-methoxy-K252d",
                                     'reason': 'Molecule too large to be a '
                                               'single amino acid'},
                                 {   'smiles': 'O=C(OC1C(O)C(OC(C1O)CO)OCC2OC(OCC(OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)C(O)C(C2O)O)CCCCCCCCCCCCCCC',
                                     'name': '[(2S)-2-hexadecanoyloxy-3-[(2S,3R,4S,5S,6R)-6-[[(2S,3R,4S,5S,6R)-4-hexadecanoyloxy-3,5-dihydroxy-6-(hydroxymethyl)tetrahydropyran-2-yl]oxymethyl]-3,4,5-trihydroxy-tetrahydropyran-2-yl]oxy-propyl] '
                                             '(9E,12E,15E)-octadeca-9,12,15-trienoate',
                                     'reason': 'Molecule too large to be a '
                                               'single amino acid'},
                                 {   'smiles': 'O=C1C2=C(O)C=C(OC)C=C2C(=O)C3=C1[C@@H]([C@@H](O)[C@]([C@@H]3O)(O)C)C',
                                     'name': 'Altersolanol G',
                                     'reason': 'Molecule too large to be a '
                                               'single amino acid'},
                                 {   'smiles': '[H][C@]1(O[C@](O)(C[C@H](O)[C@H]1NC(=O)CO)C([O-])=O)[C@H](O)[C@H](O)CO',
                                     'name': 'N-glycoloyl-alpha-neuraminate',
                                     'reason': 'Molecule too large to be a '
                                               'single amino acid'},
                                 {   'smiles': 'OC(C(O)C/C=C\\C/C=C\\C/C=C\\CC)C/C=C\\C/C=C\\CCC(O)=O',
                                     'name': '10,11-DiHDPE',
                                     'reason': 'Molecule too large to be a '
                                               'single amino acid'},
                                 {   'smiles': '[Na+].[H][C@]12SCC(C)=C(N1C(=O)[C@H]2NC(=O)[C@H](N)c1ccccc1)C([O-])=O',
                                     'name': 'cephalexin sodium',
                                     'reason': 'Molecule too large to be a '
                                               'single amino acid'}],
    'sample_false_negatives': [],
    'attempt': 1,
    'reasoning': 'Let me analyze the issues with the previous code and suggest '
                 'improvements:\n'
                 '\n'
                 '1. False Positives Analysis:\n'
                 '- Most false positives are peptides (combinations of amino '
                 'acids) rather than single amino acids\n'
                 '- The code is detecting polar groups in the entire molecule '
                 'rather than just in the side chain\n'
                 '- Some complex molecules with amino acid-like structures are '
                 'being incorrectly classified\n'
                 '\n'
                 '2. False Negatives Analysis:\n'
                 '- Tryptophan is being missed despite being polar (the indole '
                 'NH can form hydrogen bonds)\n'
                 "- The current patterns don't account for all types of "
                 'hydrogen bond donors/acceptors\n'
                 '\n'
                 '3. Key Improvements Needed:\n'
                 '- Add molecular size check to exclude peptides and larger '
                 'molecules\n'
                 '- Better isolate the side chain from the backbone for '
                 'analysis\n'
                 '- Include aromatic NH groups as polar groups\n'
                 '- Add stricter checks for being a single amino acid vs '
                 'peptide\n'
                 '\n'
                 "Here's my improved solution:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 37,
    'num_false_positives': 100,
    'num_true_negatives': 11010,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.27007299270072993,
    'recall': 1.0,
    'f1': 0.42528735632183906,
    'accuracy': 0.9910289764062079,
    'negative_predictive_value': 1.0}