"""
Classifies: CHEBI:140310 phenyl acetates
"""
"""
Classifies: phenyl acetates
Definition: An acetate ester obtained by formal condensation of the carboxy group of acetic acid 
           with the hydroxy group of any phenol.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phenyl_acetates(smiles: str):
    """
    Determines if a molecule is a phenyl acetate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) indicating if molecule is a phenyl acetate and the reason
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for phenyl acetate:
    # More specific pattern requiring phenolic oxygen and acetate group
    phenyl_acetate_pattern = Chem.MolFromSmarts("[cR1]1[cR1][cR1][cR1][cR1][cR1]1-[OX2][CX3](=[OX1])[CH3]")
    
    # Negative patterns - structures that should not be considered phenyl acetates
    negative_patterns = [
        # Complex ring systems
        Chem.MolFromSmarts("[R2]"), # Any atom in 2 or more rings
        # Specific exclusions
        Chem.MolFromSmarts("O=C1Oc2ccccc2C1"),  # coumarin
        Chem.MolFromSmarts("C1=NCc2ccccc12"),    # indole
        Chem.MolFromSmarts("O=C1C=Cc2ccccc2O1"), # benzofuranone
        Chem.MolFromSmarts("O=C1C=COc2ccccc21"), # benzopyranone
        # Exclude certain functional groups that indicate more complex structures
        Chem.MolFromSmarts("[N+](=O)[O-]"),      # nitro group
        Chem.MolFromSmarts("[#7R]"),             # nitrogen in ring
        Chem.MolFromSmarts("[OR2]"),             # oxygen in multiple rings
    ]

    # First check negative patterns
    for pattern in negative_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            return False, "Contains excluded structural features"

    # Check for phenyl acetate pattern
    matches = mol.GetSubstructMatches(phenyl_acetate_pattern)
    if not matches:
        return False, "No phenyl acetate group found"

    # Additional validation for each match
    for match in matches:
        # Get the phenyl ring atoms
        ring_atoms = match[:6]
        
        # Verify ring is a single benzene ring (not part of larger system)
        ring_info = mol.GetRingInfo()
        ring_count = 0
        for atom_idx in ring_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            # Count rings this atom is part of
            ring_count += sum(1 for ring in ring_info.AtomRings() if atom_idx in ring)
            
        # Each atom should be in exactly one ring
        if ring_count != 6:
            continue

        # Verify acetate group
        oxygen_idx = match[6]
        carbonyl_idx = match[7]
        methyl_idx = match[9]
        
        oxygen = mol.GetAtomWithIdx(oxygen_idx)
        carbonyl = mol.GetAtomWithIdx(carbonyl_idx)
        methyl = mol.GetAtomWithIdx(methyl_idx)
        
        # Verify proper connectivity
        if (len(oxygen.GetNeighbors()) != 2 or
            len(carbonyl.GetNeighbors()) != 3 or
            len(methyl.GetNeighbors()) != 1):
            continue
            
        # If we get here, we have a valid match
        return True, "Contains phenyl ring with acetate group properly attached to phenolic oxygen"

    return False, "No valid phenyl acetate group found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:140310',
                          'name': 'phenyl acetates',
                          'definition': 'An acetate ester obtained by formal '
                                        'condensation of the carboxy group of '
                                        'acetic acid with the hydroxy group of '
                                        'any phenol.',
                          'parents': ['CHEBI:22712', 'CHEBI:47622'],
                          'xrefs': ['MetaCyc:Phenyl-Acetates'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.35064935064935066 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: COc1cc(COC(C)=O)ccc1OC(C)=O NAME: '
               'vanillyl alcohol diacetate REASON: CORRECT Contains phenyl '
               'ring with acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: COC(=O)c1ccc(OC(C)=O)cc1 NAME: '
               'methyl-4-acetoxybenzoate REASON: CORRECT Contains phenyl ring '
               'with acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: COc1cc(\\C=C\\C=O)cc(OC)c1OC(C)=O NAME: '
               '4-acetoxy-3,5-dimethoxy-trans-cinnamaldehyde REASON: CORRECT '
               'Contains phenyl ring with acetate group properly attached to '
               'phenolic oxygen\n'
               ' * SMILES: CC(=O)Oc1ccccc1 NAME: phenyl acetate REASON: '
               'CORRECT Contains phenyl ring with acetate group properly '
               'attached to phenolic oxygen\n'
               ' * SMILES: COc1cc(COC(C)=O)cc(OC)c1OC(C)=O NAME: syringyl '
               'alcohol diacetate REASON: CORRECT Contains phenyl ring with '
               'acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: COc1cc(cc(OC)c1OC(C)=O)C(C)=O NAME: acetosyringone '
               'acetate REASON: CORRECT Contains phenyl ring with acetate '
               'group properly attached to phenolic oxygen\n'
               ' * SMILES: CC(=O)Oc1cccc(O)c1 NAME: resorcinol monoacetate '
               'REASON: CORRECT Contains phenyl ring with acetate group '
               'properly attached to phenolic oxygen\n'
               ' * SMILES: COC(=O)\\C=C\\c1ccc(OC(C)=O)c(OC)c1 NAME: methyl '
               '4-acetoxy-3-methoxycinnamate REASON: CORRECT Contains phenyl '
               'ring with acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: COc1cc(\\C=C\\C)ccc1OC(C)=O NAME: isoeugenol '
               'acetate REASON: CORRECT Contains phenyl ring with acetate '
               'group properly attached to phenolic oxygen\n'
               ' * SMILES: CC(=O)Oc1ccc(cc1)C(O)=O NAME: 4-acetoxy benzoic '
               'acid REASON: CORRECT Contains phenyl ring with acetate group '
               'properly attached to phenolic oxygen\n'
               ' * SMILES: CC(=O)OC1=CC=C(CO)C=C1 NAME: 4-acetoxybenzyl '
               'alcohol REASON: CORRECT Contains phenyl ring with acetate '
               'group properly attached to phenolic oxygen\n'
               ' * SMILES: CC(=O)Oc1ccc(\\C=C\\C=O)cc1 NAME: 4-acetoxy '
               'cinnamaldehyde REASON: CORRECT Contains phenyl ring with '
               'acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: COC(=O)\\C=C\\c1ccc(OC(C)=O)cc1 NAME: methyl '
               'p-coumarate acetate REASON: CORRECT Contains phenyl ring with '
               'acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: COc1cc(ccc1OC(C)=O)C(C)=O NAME: '
               '4-acetyl-2-methoxyphenyl acetate REASON: CORRECT Contains '
               'phenyl ring with acetate group properly attached to phenolic '
               'oxygen\n'
               ' * SMILES: COc1cc(C=O)cc(OC)c1OC(C)=O NAME: syringaldehyde '
               'acetate REASON: CORRECT Contains phenyl ring with acetate '
               'group properly attached to phenolic oxygen\n'
               ' * SMILES: CC(=O)OC=1C=CC(=CC1OC)C(=O)OC NAME: methyl '
               '4-acetoxy-3-methoxybenzoate REASON: CORRECT Contains phenyl '
               'ring with acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: COc1cc(\\C=C\\C=O)ccc1OC(C)=O NAME: '
               '4-acetoxy-3-methoxy-trans-cinnamaldehyde REASON: CORRECT '
               'Contains phenyl ring with acetate group properly attached to '
               'phenolic oxygen\n'
               ' * SMILES: CC(=O)Oc1ccc(O)cc1 NAME: 4-hydroxyphenyl acetate '
               'REASON: CORRECT Contains phenyl ring with acetate group '
               'properly attached to phenolic oxygen\n'
               ' * SMILES: CC(=O)Nc1cc(C)ccc1OC(C)=O NAME: '
               '2-(acetylamino)-4-methylphenyl acetate REASON: CORRECT '
               'Contains phenyl ring with acetate group properly attached to '
               'phenolic oxygen\n'
               ' * SMILES: CC(=O)OC\\C=C\\c1ccc(OC(C)=O)cc1 NAME: 4-hydroxy '
               'cinnamyl alcohol diacetate REASON: CORRECT Contains phenyl '
               'ring with acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: CC(=O)Oc1ccc(C=O)cc1 NAME: 4-acetoxy benzaldehyde '
               'REASON: CORRECT Contains phenyl ring with acetate group '
               'properly attached to phenolic oxygen\n'
               ' * SMILES: COc1cc(\\C=C\\COC(C)=O)ccc1OC(C)=O NAME: '
               'trans-coniferyl alcohol diacetate REASON: CORRECT Contains '
               'phenyl ring with acetate group properly attached to phenolic '
               'oxygen\n'
               ' * SMILES: CC(=O)Oc1ccc(cc1)C(C)=O NAME: 4-acetoxy '
               'acetophenone REASON: CORRECT Contains phenyl ring with acetate '
               'group properly attached to phenolic oxygen\n'
               ' * SMILES: COc1cc(\\C=C\\C(O)=O)ccc1OC(C)=O NAME: '
               'O-acetylferulic acid REASON: CORRECT Contains phenyl ring with '
               'acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: COc1cc(cc(OC)c1OC(C)=O)C(O)=O NAME: syringic acid '
               'acetate REASON: CORRECT Contains phenyl ring with acetate '
               'group properly attached to phenolic oxygen\n'
               ' * SMILES: CC(=O)Oc1ccc(cc1)[N+]([O-])=O NAME: 4-nitrophenyl '
               'acetate REASON: CORRECT Contains phenyl ring with acetate '
               'group properly attached to phenolic oxygen\n'
               ' * SMILES: COc1cc(C=O)ccc1OC(C)=O NAME: vanillin acetate '
               'REASON: CORRECT Contains phenyl ring with acetate group '
               'properly attached to phenolic oxygen\n'
               'False positives: SMILES: '
               'C(OC1=CC=C(C=C1)NC(C)=O)(=O)C2=C(OC(C)=O)C=CC=C2 NAME: '
               'benorilate REASON: WRONGLY CLASSIFIED Contains phenyl ring '
               'with acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: CC(=O)OC1=CC=CC=C1C(=O)OC2=CC=CC=C2C(=O)O NAME: '
               '2-[(2-acetyloxyphenyl)-oxomethoxy]benzoic acid REASON: WRONGLY '
               'CLASSIFIED Contains phenyl ring with acetate group properly '
               'attached to phenolic oxygen\n'
               ' * SMILES: '
               'O=C(OC1=C(OC(=O)C)C(=C(O)C(=C1C2=CC=C(O)C=C2)O)C3=CC=C(O)C=C3)C4=CC=CC=C4 '
               'NAME: Concrescenin A REASON: WRONGLY CLASSIFIED Contains '
               'phenyl ring with acetate group properly attached to phenolic '
               'oxygen\n'
               ' * SMILES: '
               'O=C(OC1=C(O)C(=C(OC(=O)C)C(=C1C2=CC=C(O)C=C2)O)C3=CC=C(O)C=C3)C '
               'NAME: Kynapcin-12 REASON: WRONGLY CLASSIFIED Contains phenyl '
               'ring with acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: O(C1=C(OC)C=C(OC(=O)C)C=C1OC)C NAME: '
               '3,4,5-Trimethoxyphenyl acetate REASON: WRONGLY CLASSIFIED '
               'Contains phenyl ring with acetate group properly attached to '
               'phenolic oxygen\n'
               ' * SMILES: CC(=O)OC1=C(C=C(C=C1)C(=O)OC2CC3CCC2(C3(C)C)C)OC '
               'NAME: 4-acetyloxy-3-methoxybenzoic acid '
               '(4,7,7-trimethyl-3-bicyclo[2.2.1]heptanyl) ester REASON: '
               'WRONGLY CLASSIFIED Contains phenyl ring with acetate group '
               'properly attached to phenolic oxygen\n'
               ' * SMILES: '
               'O=C(OC1=C(OC(=O)C)C(=C(O)C(=C1C2=CC(O)=C(O)C=C2)O)C3=CC=C(OC(=O)C)C=C3)C '
               'NAME: 2′,3′,4′′-triacetoxy-3,4,5′,6′-tetrahydroxy-p-terphenyl '
               'REASON: WRONGLY CLASSIFIED Contains phenyl ring with acetate '
               'group properly attached to phenolic oxygen\n'
               ' * SMILES: '
               'O=C1OC2=C(C=CC(=C2)C3=C(OC(=O)C)C(OC(=O)C)=C(C4=CC=C(OC(=O)C)C=C4)C(=C3O)O)OC1(N(O)C(=O)/C(=N/O)/[C@H](CC)C)[C@H](CC)C '
               'NAME: Sarcodonin REASON: WRONGLY CLASSIFIED Contains phenyl '
               'ring with acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: '
               'CC(=O)OC1=CC(=NC2=CC=CC=C21)C=CC3=CC(=CC(=C3OC(=O)C)[N+](=O)[O-])OC '
               'NAME: acetic acid '
               '[2-[2-(2-acetyloxy-5-methoxy-3-nitrophenyl)ethenyl]-4-quinolinyl] '
               'ester REASON: WRONGLY CLASSIFIED Contains phenyl ring with '
               'acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: CCC(=O)C(CC)(c1ccc(OC(C)=O)cc1)c1ccc(OC(C)=O)cc1 '
               'NAME: 4,4-Bis[4-(acetyloxy)phenyl]3-hexanone REASON: WRONGLY '
               'CLASSIFIED Contains phenyl ring with acetate group properly '
               'attached to phenolic oxygen\n'
               ' * SMILES: O(C=1C(=CC=CC1)C)C(=O)C NAME: o-Tolyl acetate '
               'REASON: WRONGLY CLASSIFIED Contains phenyl ring with acetate '
               'group properly attached to phenolic oxygen\n'
               ' * SMILES: '
               'O=C1C(C#CC(=C)CC/C=C(/C(=O)O)\\C)=C[C@@H](OC(=O)C)[C@@H]([C@H]1OC(=O)/C(=C/CCC(C#CC2=C(O)C(O)=CC(=C2)OC(=O)C)=C)/C)O '
               'NAME: Tricholomenyn E REASON: WRONGLY CLASSIFIED Contains '
               'phenyl ring with acetate group properly attached to phenolic '
               'oxygen\n'
               ' * SMILES: CC(=O)OC1=CC=CC=C1C(=O)NC2=CC=CC=C2 NAME: acetic '
               'acid [2-[anilino(oxo)methyl]phenyl] ester REASON: WRONGLY '
               'CLASSIFIED Contains phenyl ring with acetate group properly '
               'attached to phenolic oxygen\n'
               ' * SMILES: '
               'O=C(OC1=C(O)C(=C(OC(=O)C)C(=C1C2=CC=C(O)C=C2)O)C3=CC=C(O)C=C3)CCC '
               'NAME: Curtisian O REASON: WRONGLY CLASSIFIED Contains phenyl '
               'ring with acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: '
               'CC(O)CC(=O)Oc1c(OC(=O)CCc2ccccc2)c(-c2ccc(O)cc2)c(OC(C)=O)c(OC(C)=O)c1-c1ccc(O)cc1 '
               'NAME: curtisian D REASON: WRONGLY CLASSIFIED Contains phenyl '
               'ring with acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: CC(=O)N1C(CC(=N1)C2=CC=CC=C2OC(=O)C)C3=CC=C(C=C3)OC '
               'NAME: acetic acid '
               '[2-[2-acetyl-3-(4-methoxyphenyl)-3,4-dihydropyrazol-5-yl]phenyl] '
               'ester REASON: WRONGLY CLASSIFIED Contains phenyl ring with '
               'acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: '
               'O=C(OC1=C(O)C(=C(OC(=O)C)C(=C1C2=CC=C(O)C=C2)O)C3=CC=C(O)C=C3)CCC4=CC=CC=C4 '
               'NAME: Curtisian P REASON: WRONGLY CLASSIFIED Contains phenyl '
               'ring with acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: CCn1c(c(C)c2cc(OC(C)=O)ccc12)-c1ccc(OC(C)=O)cc1 '
               'NAME: Zindoxifene REASON: WRONGLY CLASSIFIED Contains phenyl '
               'ring with acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: '
               'O=C(OC1=C(OC(=O)C)C(=C(O)C(=C1C2=CC=C(O)C=C2)O)C3=CC=C(O)C=C3)CC4=CC=CC=C4 '
               'NAME: Terrestrin D REASON: WRONGLY CLASSIFIED Contains phenyl '
               'ring with acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: CC(=O)OC1=C(C=C(C=C1)C(=O)NC2=CC=C(C=C2)O)OC NAME: '
               'acetic acid [4-[(4-hydroxyanilino)-oxomethyl]-2-methoxyphenyl] '
               'ester REASON: WRONGLY CLASSIFIED Contains phenyl ring with '
               'acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: O(C1=C(OC)C=C(OC)C=C1OC)C(=O)C NAME: '
               '2,4,6-Trimethoxyphenyl acetate REASON: WRONGLY CLASSIFIED '
               'Contains phenyl ring with acetate group properly attached to '
               'phenolic oxygen\n'
               ' * SMILES: C1(=C(C=C(OC(C)=O)C(=C1)C)C(C)C)OCCN(C)C.Cl NAME: '
               'Moxisylyte hydrochloride REASON: WRONGLY CLASSIFIED Contains '
               'phenyl ring with acetate group properly attached to phenolic '
               'oxygen\n'
               ' * SMILES: S(=O)(=O)(OC1=C(C(OC(=O)C)=CC(=C1)CCCCC)CCCC)O '
               'NAME: Stemphol B REASON: WRONGLY CLASSIFIED Contains phenyl '
               'ring with acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: OC(C1=CC=C(OC(=O)C)C=C1)C=C NAME: '
               "1'-Hydroxychavicol acetate REASON: WRONGLY CLASSIFIED Contains "
               'phenyl ring with acetate group properly attached to phenolic '
               'oxygen\n'
               ' * SMILES: O(C=1C(OC)=CC=CC1)C(=O)C NAME: Guaicyl acetate '
               'REASON: WRONGLY CLASSIFIED Contains phenyl ring with acetate '
               'group properly attached to phenolic oxygen\n'
               ' * SMILES: O(C=1C=C(CC=C)C=CC1OC(=O)C)C(=O)C NAME: '
               '4-Allyl-1,2-diacetoxybenzene REASON: WRONGLY CLASSIFIED '
               'Contains phenyl ring with acetate group properly attached to '
               'phenolic oxygen\n'
               ' * SMILES: O(C=1C(=C(CCCCCCCCCCCCCCC)C=CC1)C(O)=O)C(=O)C NAME: '
               '2-Acetoxy-6-pentadecylbenzoic acid REASON: WRONGLY CLASSIFIED '
               'Contains phenyl ring with acetate group properly attached to '
               'phenolic oxygen\n'
               ' * SMILES: O(C1=C(C=CC=C1C(O)=O)C)C(=O)C NAME: CRESOPYRINE '
               'REASON: WRONGLY CLASSIFIED Contains phenyl ring with acetate '
               'group properly attached to phenolic oxygen\n'
               ' * SMILES: CCCCCC\\C=C/CCCCCCCc1cc(O)cc(O)c1OC(C)=O NAME: '
               'ardisiphenol B REASON: WRONGLY CLASSIFIED Contains phenyl ring '
               'with acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: CC(=O)OC1=CC(=C(C=C1)NC(=O)C2CCCCC2)C(=O)O NAME: '
               '5-acetyloxy-2-[[cyclohexyl(oxo)methyl]amino]benzoic acid '
               'REASON: WRONGLY CLASSIFIED Contains phenyl ring with acetate '
               'group properly attached to phenolic oxygen\n'
               ' * SMILES: O(C1=C(C(=C(O)C=C1C)C)C)C(=O)C NAME: '
               '4-Hydroxy-2,3,6-trimethylphenyl acetate REASON: WRONGLY '
               'CLASSIFIED Contains phenyl ring with acetate group properly '
               'attached to phenolic oxygen\n'
               ' * SMILES: '
               'O1[C@@]([C@H](O)[C@@H](O)[C@H](O)[C@H]1COC(=O)C)(C=2C=3OC(CC(=O)C3C=CC2O)C4=CC=C(OC(=O)C)C=C4)[H] '
               "NAME: 4',6''-Diacetyl puerarin REASON: WRONGLY CLASSIFIED "
               'Contains phenyl ring with acetate group properly attached to '
               'phenolic oxygen\n'
               ' * SMILES: '
               'CC(=O)Oc1c(OC(C)=O)c(-c2ccc(O)cc2)c(OC(=O)c2ccccc2)c(OC(C)=O)c1-c1ccc(O)cc1 '
               'NAME: curtisian A REASON: WRONGLY CLASSIFIED Contains phenyl '
               'ring with acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: '
               'CC(=O)Oc1c(OC(C)=O)c(-c2ccc(O)cc2)c(OC(=O)CCc2ccccc2)c(OC(C)=O)c1-c1ccc(O)cc1 '
               'NAME: curtisian B REASON: WRONGLY CLASSIFIED Contains phenyl '
               'ring with acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: O=C1OC(C2=CC(=C(OC(=O)C)C(=C2)C)CC)=CC=C1C NAME: '
               'Penipyrol B REASON: WRONGLY CLASSIFIED Contains phenyl ring '
               'with acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: '
               'CCOC1=C(C=CC(=C1[N+](=O)[O-])C=CC2=NC3=C(C=CC=C3OC(=O)C)C=C2)OC(=O)C '
               'NAME: acetic acid '
               '[4-[2-(8-acetyloxy-2-quinolinyl)ethenyl]-2-ethoxy-3-nitrophenyl] '
               'ester REASON: WRONGLY CLASSIFIED Contains phenyl ring with '
               'acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: '
               'C=12C(=CC=CC1)N=C(S2)C=3C=C(C=CC3OC(C=4C(=CC=CC4)OC(C)=O)=O)F '
               'NAME: aspirin-based probe AP REASON: WRONGLY CLASSIFIED '
               'Contains phenyl ring with acetate group properly attached to '
               'phenolic oxygen\n'
               ' * SMILES: '
               'O=C(OC1=CC(=C([C@@H]([C@@H](O)[C@@H](C(=O)C)C)C)C=C1)/C=C/CO)C '
               'NAME: Benwamycin C REASON: WRONGLY CLASSIFIED Contains phenyl '
               'ring with acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: '
               'CC[C@H](C)C1=C(O)[N+]2([O-])Oc3ccc(cc3O[C@]2([C@@H](C)CC)C(=O)N1OC)-c1c(O)c(O)c(-c2ccc(O)cc2)c(OC(C)=O)c1OC(C)=O '
               'NAME: Sarcodonin-delta REASON: WRONGLY CLASSIFIED Contains '
               'phenyl ring with acetate group properly attached to phenolic '
               'oxygen\n'
               ' * SMILES: '
               'O=C(OC1=C(OC(=O)C2=CC=CC=C2)C(=C(OC(=O)C)C(=C1C3=CC=C(O)C=C3)OC(=O)C4=CC=CC=C4)C5=CC=C(O)C=C5)C6=CC=CC=C6 '
               'NAME: Thelephantin L REASON: WRONGLY CLASSIFIED Contains '
               'phenyl ring with acetate group properly attached to phenolic '
               'oxygen\n'
               ' * SMILES: O(C(C)=O)C1=C(C(OC2=CC=CC=C2)=O)C=CC=C1 NAME: '
               'phenyl acetylsalicylate REASON: WRONGLY CLASSIFIED Contains '
               'phenyl ring with acetate group properly attached to phenolic '
               'oxygen\n'
               ' * SMILES: O(C1=CC=C(CCC(=O)C)C=C1)C(=O)C NAME: Cuelure '
               'REASON: WRONGLY CLASSIFIED Contains phenyl ring with acetate '
               'group properly attached to phenolic oxygen\n'
               ' * SMILES: O=C(CCCCCCCCC)CCC1=CC(OC)=C(OC(=O)C)C=C1 NAME: '
               '[8]-Paradyl acetate REASON: WRONGLY CLASSIFIED Contains phenyl '
               'ring with acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: '
               'C(C1=CC=C(C=C1)OC(C)=O)(C2=CC=C(C=C2)OC(C)=O)=C3CCCCC3 NAME: '
               'Cyclofenil REASON: WRONGLY CLASSIFIED Contains phenyl ring '
               'with acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: '
               'O(C=1C(C/C=C(\\CC\\C=C(\\CC\\C=C(\\CCC=C(C)C)/C)/C)/C)=C(O)C=CC1O)C(=O)C '
               'NAME: 2-Acetoxy-3-geranylgeranyl-1,4-dihydroxybenzene REASON: '
               'WRONGLY CLASSIFIED Contains phenyl ring with acetate group '
               'properly attached to phenolic oxygen\n'
               ' * SMILES: '
               'CCOC1=CC(=CC(=C1OC(=O)C)[N+](=O)[O-])C=C2C(=O)SC(=N2)SCC=C '
               'NAME: acetic acid '
               '[2-ethoxy-6-nitro-4-[[5-oxo-2-(prop-2-enylthio)-4-thiazolylidene]methyl]phenyl] '
               'ester REASON: WRONGLY CLASSIFIED Contains phenyl ring with '
               'acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: CC(=O)OC1=CC=CC=C1C(=O)NC2=NC=C(S2)[N+](=O)[O-] '
               'NAME: acetic acid '
               '[2-[[(5-nitro-2-thiazolyl)amino]-oxomethyl]phenyl] ester '
               'REASON: WRONGLY CLASSIFIED Contains phenyl ring with acetate '
               'group properly attached to phenolic oxygen\n'
               ' * SMILES: '
               'O=C1O[C@](C(=O)OC)(CC2=CC3=C(OC(C)(C)C(C3)OC(=O)C)C=C2)C(=C1O)C4=CC=C(OC(=O)C)C=C4 '
               "NAME: 4',8''-diacetoxy butyrolactone VI REASON: WRONGLY "
               'CLASSIFIED Contains phenyl ring with acetate group properly '
               'attached to phenolic oxygen\n'
               ' * SMILES: CC(=O)OC1=CC=CC=C1C(=O)NC2=CC(=O)C=CC2(OC)OC NAME: '
               'acetic acid '
               '[2-[[(6,6-dimethoxy-3-oxo-1-cyclohexa-1,4-dienyl)amino]-oxomethyl]phenyl] '
               'ester REASON: WRONGLY CLASSIFIED Contains phenyl ring with '
               'acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: '
               'O(C1=C(C(OC(=O)C)=C(OC(=O)C)C(=C1OC(=O)C)C2=CC=C(OC(=O)C)C=C2)C3=CC(OC(=O)C)=C(OC(=O)C)C=C3)C(=O)C '
               'NAME: Protoleucomelone REASON: WRONGLY CLASSIFIED Contains '
               'phenyl ring with acetate group properly attached to phenolic '
               'oxygen\n'
               ' * SMILES: '
               'CC(=O)OC1=CC=C(C=C1)C2=CC=C(C=C2)C(=O)OCC(=O)N(C)CC(=O)NC3=CC=CC=C3Cl '
               'NAME: 4-(4-acetyloxyphenyl)benzoic acid '
               '[2-[[2-(2-chloroanilino)-2-oxoethyl]-methylamino]-2-oxoethyl] '
               'ester REASON: WRONGLY CLASSIFIED Contains phenyl ring with '
               'acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: '
               'O=C(OC1=C(C(O)=C2OC3=C(C2=C1OC(=O)C)C=C(O)C(=C3)O)C4=CC=C(OC(=O)C)C=C4)C '
               'NAME: Sarcodan REASON: WRONGLY CLASSIFIED Contains phenyl ring '
               'with acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: CC(C)c1ccc(C)cc1OC(C)=O NAME: Thymyl acetate '
               'REASON: WRONGLY CLASSIFIED Contains phenyl ring with acetate '
               'group properly attached to phenolic oxygen\n'
               ' * SMILES: CC(=O)N1CCC(=CC1)C2=CC(=C(C=C2)OC(=O)C)OC NAME: '
               'acetic acid '
               '[4-(1-acetyl-3,6-dihydro-2H-pyridin-4-yl)-2-methoxyphenyl] '
               'ester REASON: WRONGLY CLASSIFIED Contains phenyl ring with '
               'acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: C(OC1=C(OC)C=CC=C1)(=O)C2=C(OC(C)=O)C=CC=C2 NAME: '
               'guacetisal REASON: WRONGLY CLASSIFIED Contains phenyl ring '
               'with acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: CC1=CC(=C(C=C1OC(=O)C)C(C)C)OCCN(C)C NAME: acetic '
               'acid '
               '[4-[2-(dimethylamino)ethoxy]-2-methyl-5-propan-2-ylphenyl] '
               'ester REASON: WRONGLY CLASSIFIED Contains phenyl ring with '
               'acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: '
               'O=C(OC1=C(OC(=O)C)C(=C(O)C(=C1C2=CC(O)=C(O)C=C2)O)C3=CC=C(O)C=C3)C '
               "NAME: 2',3'-diacetoxy-3,4,5',6',4''-pentahydroxy-p-terphenyl "
               'REASON: WRONGLY CLASSIFIED Contains phenyl ring with acetate '
               'group properly attached to phenolic oxygen\n'
               ' * SMILES: CC1=CC=C(C=C1)S(=O)(=O)N(C2=CC=CC=C2OC(=O)C)C(=O)C '
               'NAME: acetic acid '
               '[2-[acetyl-(4-methylphenyl)sulfonylamino]phenyl] ester REASON: '
               'WRONGLY CLASSIFIED Contains phenyl ring with acetate group '
               'properly attached to phenolic oxygen\n'
               ' * SMILES: '
               'O1[C@@H]([C@H](OC(=O)C2=CC(OC(=O)C)=C(OC(=O)C)C(OC(=O)C)=C2)CC3=C1C=C(OC(=O)C)C=C3OC(=O)C)C4=CC(OC(=O)C)=C(OC(=O)C)C(OC(=O)C)=C4 '
               'NAME: EGCG Octaacetate REASON: WRONGLY CLASSIFIED Contains '
               'phenyl ring with acetate group properly attached to phenolic '
               'oxygen\n'
               ' * SMILES: '
               'CCOC(=O)C1=C(N=C2N(C1C3=CC(=C(C=C3)OC(=O)C)OC)C(=O)C(S2)C)C '
               'NAME: '
               '5-(4-acetyloxy-3-methoxyphenyl)-2,7-dimethyl-3-oxo-5H-thiazolo[3,2-a]pyrimidine-6-carboxylic '
               'acid ethyl ester REASON: WRONGLY CLASSIFIED Contains phenyl '
               'ring with acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: '
               'CC(=O)OC1=CC=C(C=C1)C(=O)NC2=CC=CC(=C2)C3=NN=C(O3)C4=CC=CO4 '
               'NAME: acetic acid '
               '[4-[[3-[5-(2-furanyl)-1,3,4-oxadiazol-2-yl]anilino]-oxomethyl]phenyl] '
               'ester REASON: WRONGLY CLASSIFIED Contains phenyl ring with '
               'acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: '
               'O=C(OC1=C(O)C(=C(OC(=O)C)C(=C1C2=CC=C(O)C=C2)O)C3=CC=C(O)C=C3)C[C@@H](OC(=O)C)C '
               'NAME: Curtisian N REASON: WRONGLY CLASSIFIED Contains phenyl '
               'ring with acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: CCCCC\\C=C/C\\C=C/CCCCCCCc1cc(O)cc(O)c1OC(C)=O '
               'NAME: ardisiphenol C REASON: WRONGLY CLASSIFIED Contains '
               'phenyl ring with acetate group properly attached to phenolic '
               'oxygen\n'
               ' * SMILES: CC(=O)O[C@@H](C=C)c1ccc(OC(C)=O)cc1 NAME: '
               "1'-acetoxychavicol acetate REASON: WRONGLY CLASSIFIED Contains "
               'phenyl ring with acetate group properly attached to phenolic '
               'oxygen\n'
               ' * SMILES: '
               'COc1cc(cc(OC)c1OC(C)=O)C1OCC2C1COC2c1cc(OC)c(OC(C)=O)c(OC)c1 '
               'NAME: syringylresinol diacetate REASON: WRONGLY CLASSIFIED '
               'Contains phenyl ring with acetate group properly attached to '
               'phenolic oxygen\n'
               ' * SMILES: '
               'O(C1=C(C\\C=C(\\CC\\C=C(\\CC\\C=C(\\CCC=C(C)C)/C)/C)/C)C(O)=C(O)C=C1)C(=O)C '
               'NAME: Suillin REASON: WRONGLY CLASSIFIED Contains phenyl ring '
               'with acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: '
               'O=C(OC1=C(OC(=O)CCC2=CC=CC=C2)C(=C(OC(=O)C[C@@H](O)C)C(=C1C3=CC=C(O)C=C3)OC(=O)C)C4=CC=C(O)C=C4)C[C@@H](O)C '
               'NAME: Curtisian E REASON: WRONGLY CLASSIFIED Contains phenyl '
               'ring with acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: COc1cc(CC=C)ccc1OC(C)=O NAME: Acetyleugenol REASON: '
               'WRONGLY CLASSIFIED Contains phenyl ring with acetate group '
               'properly attached to phenolic oxygen\n'
               ' * SMILES: CCCCC#CCSC1=CC=CC=C1OC(C)=O NAME: '
               '(2-hept-2-ynylsulfanylphenyl) acetate REASON: WRONGLY '
               'CLASSIFIED Contains phenyl ring with acetate group properly '
               'attached to phenolic oxygen\n'
               ' * SMILES: '
               'O=C(OC1=C(OC(=O)C)C(=C(OC(=O)C2=CC=CC=C2)C(=C1C3=CC=C(O)C=C3)OC(=O)C)C4=CC=C(O)C=C4)C5=CC=CC=C5 '
               'NAME: Concrescenin B REASON: WRONGLY CLASSIFIED Contains '
               'phenyl ring with acetate group properly attached to phenolic '
               'oxygen\n'
               ' * SMILES: '
               'CC(=O)OC1=CC=CC=C1C2=NN(C(C2)C3=CC=CC=C3)C(=O)C=CC4=CC=CC=C4 '
               'NAME: acetic acid '
               '[2-[2-(1-oxo-3-phenylprop-2-enyl)-3-phenyl-3,4-dihydropyrazol-5-yl]phenyl] '
               'ester REASON: WRONGLY CLASSIFIED Contains phenyl ring with '
               'acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: '
               'CC(=O)OC1=C(C=C(C=C1)C(=O)NC2=C(C=C(C=C2)O)C(=O)O)OC NAME: '
               '2-[[(4-acetyloxy-3-methoxyphenyl)-oxomethyl]amino]-5-hydroxybenzoic '
               'acid REASON: WRONGLY CLASSIFIED Contains phenyl ring with '
               'acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: CC(=O)OC1=CC=C(C=C1)NC(=O)CC2=CC=CC=C2 NAME: acetic '
               'acid [4-[(1-oxo-2-phenylethyl)amino]phenyl] ester REASON: '
               'WRONGLY CLASSIFIED Contains phenyl ring with acetate group '
               'properly attached to phenolic oxygen\n'
               ' * SMILES: '
               '[OH-].[OH-].[Mg++].[Al+3].[O-]C([O-])=O.NCC([O-])=O.CC(=O)Oc1ccccc1C(O)=O '
               'NAME: Bufferin REASON: WRONGLY CLASSIFIED Contains phenyl ring '
               'with acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: CC(=O)Oc1ccc(CCC(O)=O)cc1 NAME: '
               "3-(4'-acetoxyphenyl)propionic acid REASON: WRONGLY CLASSIFIED "
               'Contains phenyl ring with acetate group properly attached to '
               'phenolic oxygen\n'
               ' * SMILES: CC(=O)NC1=NC(=NC(=N1)C2=CC=CC=C2OC(=O)C)C3=CC=CC=C3 '
               'NAME: acetic acid '
               '[2-(4-acetamido-6-phenyl-1,3,5-triazin-2-yl)phenyl] ester '
               'REASON: WRONGLY CLASSIFIED Contains phenyl ring with acetate '
               'group properly attached to phenolic oxygen\n'
               ' * SMILES: '
               'O=C(OC1=C(OC(=O)C[C@H](O)C)C(=C(OC(=O)C)C(=C1C2=CC=C(O)C=C2)OC(=O)CCC3=CC=CC=C3)C4=CC=C(O)C=C4)C[C@H](O)C '
               'NAME: Curtisian V REASON: WRONGLY CLASSIFIED Contains phenyl '
               'ring with acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: CC(=O)OC1=CC=C(C=C1)C(=O)NCCC2=CC=CC=C2 NAME: '
               'acetic acid [4-[oxo-(2-phenylethylamino)methyl]phenyl] ester '
               'REASON: WRONGLY CLASSIFIED Contains phenyl ring with acetate '
               'group properly attached to phenolic oxygen\n'
               ' * SMILES: '
               'CC(=O)Oc1c(cc(c(C)c1C(C)(C)C)[N+]([O-])=O)[N+]([O-])=O NAME: '
               'Medinoterb acetate REASON: WRONGLY CLASSIFIED Contains phenyl '
               'ring with acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: CC(=O)OC1=CC=CC=C1C(=O)OC2=CC=CC(=C2)CO[N+](=O)[O-] '
               'NAME: 2-acetyloxybenzoic acid [3-(nitrooxymethyl)phenyl] ester '
               'REASON: WRONGLY CLASSIFIED Contains phenyl ring with acetate '
               'group properly attached to phenolic oxygen\n'
               ' * SMILES: '
               'CC(O)CC(=O)Oc1c(OC(=O)CC(C)OC(C)=O)c(-c2ccc(O)cc2)c(OC(C)=O)c(OC(=O)CC(C)OC(C)=O)c1-c1ccc(O)cc1 '
               'NAME: curtisian C REASON: WRONGLY CLASSIFIED Contains phenyl '
               'ring with acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: '
               'C(N(C)C1CCCCC1)C2=C(NC(=O)C3=CC(OC)=C(C=C3)OC(C)=O)C(Br)=CC(=C2)Br '
               'NAME: brovanexine REASON: WRONGLY CLASSIFIED Contains phenyl '
               'ring with acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: CC(=O)OC1=CC=CC(=C1OC(=O)C)C(=O)O NAME: '
               '2,3-diacetyloxybenzoic acid REASON: WRONGLY CLASSIFIED '
               'Contains phenyl ring with acetate group properly attached to '
               'phenolic oxygen\n'
               ' * SMILES: '
               'O(C=1C(CCCCCCCCCCC(=O)C2=C(O)C=C(O)C=C2O)=C(O)C=CC1)C(=O)C '
               'NAME: Ardisinone A REASON: WRONGLY CLASSIFIED Contains phenyl '
               'ring with acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: '
               'CCOC1=C(C=CC(=C1)C=C2C(=O)N3C(C(=C(N=C3S2)C)C(=O)OCC)C=CC4=CC=CC=C4)OC(=O)C '
               'NAME: '
               '2-[(4-acetyloxy-3-ethoxyphenyl)methylidene]-7-methyl-3-oxo-5-(2-phenylethenyl)-5H-thiazolo[3,2-a]pyrimidine-6-carboxylic '
               'acid ethyl ester REASON: WRONGLY CLASSIFIED Contains phenyl '
               'ring with acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: CC(=O)OC1=C(C=C(C=C1)C(=O)NC2=CC=C(C=C2)C(=O)O)OC '
               'NAME: '
               '4-[[(4-acetyloxy-3-methoxyphenyl)-oxomethyl]amino]benzoic acid '
               'REASON: WRONGLY CLASSIFIED Contains phenyl ring with acetate '
               'group properly attached to phenolic oxygen\n'
               ' * SMILES: O=C(OC1=C(OC2=C(C=C(OC(=O)C)C=C2C)C)C(=CC(=C1)C)C)C '
               'NAME: Aspergillether A REASON: WRONGLY CLASSIFIED Contains '
               'phenyl ring with acetate group properly attached to phenolic '
               'oxygen\n'
               ' * SMILES: '
               'CC(=O)OC1=CC=C(C=C1)C=CC2=CC(=CC(=C2)OC(=O)C)OC(=O)C NAME: '
               'acetic acid [4-[2-(3,5-diacetyloxyphenyl)ethenyl]phenyl] ester '
               'REASON: WRONGLY CLASSIFIED Contains phenyl ring with acetate '
               'group properly attached to phenolic oxygen\n'
               ' * SMILES: '
               'O(C=1C(C(N2CCCCC2)=CC=C(C3=CC=CC=C3)C4=CC=CC=C4)=CC=CC1)C(=O)C '
               'NAME: 2-(4,4-diphenyl-1-piperidinobuta-1,3-dienyl)phenyl '
               'acetate REASON: WRONGLY CLASSIFIED Contains phenyl ring with '
               'acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: CCC(C)c1cc(cc(c1OC(C)=O)[N+]([O-])=O)[N+]([O-])=O '
               'NAME: Dinoseb acetate REASON: WRONGLY CLASSIFIED Contains '
               'phenyl ring with acetate group properly attached to phenolic '
               'oxygen\n'
               ' * SMILES: '
               'CC(=O)OC1=CC=C(C=C1)C(C2=CC=C(C=C2)OC(=O)C)C3=CC=CC=N3 NAME: '
               'Bisacodyl REASON: WRONGLY CLASSIFIED Contains phenyl ring with '
               'acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: '
               'O1C([C@@H](OC(=O)C)[C@H](OC(=O)C)C(OC(=O)C)[C@@H]1OC=2C=C3OC(CC(=O)C3=C(O)C2)C4=CC(OC(=O)C)=C(OC(=O)C)C=C4)COC(=O)C '
               'NAME: Hexaacetylpyracanthoside REASON: WRONGLY CLASSIFIED '
               'Contains phenyl ring with acetate group properly attached to '
               'phenolic oxygen\n'
               ' * SMILES: CCC1=C(c2ccc(OC(C)=O)cc2)C(C)(C)Oc2cc(OC(C)=O)ccc12 '
               'NAME: '
               '2,2-dimethyl-3-[4-(acetyloxy)phenyl]-4-ethyl-2H-1-benzopyran-7-ol '
               'acetate REASON: WRONGLY CLASSIFIED Contains phenyl ring with '
               'acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: O(C1=CC=C(CC=C)C=C1)C(=O)C NAME: Chavicol acetate '
               'REASON: WRONGLY CLASSIFIED Contains phenyl ring with acetate '
               'group properly attached to phenolic oxygen\n'
               ' * SMILES: '
               'O=C1C(C2=C(N1)C=CC=C2)(C3=CC=C(C=C3)OC(C)=O)C4=CC=C(C=C4)OC(C)=O '
               'NAME: oxyphenisatine acetate REASON: WRONGLY CLASSIFIED '
               'Contains phenyl ring with acetate group properly attached to '
               'phenolic oxygen\n'
               ' * SMILES: CC(C)NCC(O)COc1cc(C)c(OC(C)=O)c(C)c1C NAME: '
               'metipranolol REASON: WRONGLY CLASSIFIED Contains phenyl ring '
               'with acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: '
               'O1C2=C(C=3C1=CC(O)=C(O)C3)C(OC(=O)C)=C(OC(=O)C)C(=C2OC(=O)C)C4=CC=C(OC(=O)C)C=C4 '
               'NAME: BL III REASON: WRONGLY CLASSIFIED Contains phenyl ring '
               'with acetate group properly attached to phenolic oxygen\n'
               ' * SMILES: O=C(OC1=CC(OC2=CC(OC(=O)C)=CC(=C2)C)=CC(=C1)C)C '
               'NAME: Trichoderolide A REASON: WRONGLY CLASSIFIED Contains '
               'phenyl ring with acetate group properly attached to phenolic '
               'oxygen\n'
               ' * SMILES: '
               'O=C1C(OC2=C(N1)C=CC=C2)(C3=CC=C(C=C3)OC(C)=O)C4=CC=C(C=C4)OC(C)=O '
               'NAME: Bisoxatin acetate REASON: WRONGLY CLASSIFIED Contains '
               'phenyl ring with acetate group properly attached to phenolic '
               'oxygen\n'
               ' * SMILES: '
               'O=C(OC1=CC(O)=C(C(=O)OC2=CC(O)=C(C(=O)O)C(=C2)C)C(=C1)C)C3=C(OC)C=C(O)C(=C3C)OC(=O)C '
               'NAME: Deliseic acid REASON: WRONGLY CLASSIFIED Contains phenyl '
               'ring with acetate group properly attached to phenolic oxygen\n'
               'False negatives: NONE\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'O=C1O[C@H]2[C@H](O[C@H](C2)C[C@H](O)C)C=3C1=C(O)C(OC)=C(OC)C3',
                                     'name': '(12R)-12-hydroxymonocerin',
                                     'reason': 'Contains excluded structural '
                                               'features'},
                                 {   'smiles': 'P(OC[C@H](OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCCCC/C=C\\CCCCCCCC)(OCCN)(O)=O',
                                     'name': 'PE(20:1(11Z)/16:0)',
                                     'reason': 'No phenyl acetate group found'},
                                 {   'smiles': 'O=C(C1=C(N)C2=C(OC(C)(C)C=C2)C=C1)C[C@@H](NC(=O)C)CO',
                                     'name': 'Fusarochromene',
                                     'reason': 'Contains excluded structural '
                                               'features'},
                                 {   'smiles': 'COc1cc(\\C=C\\C(=O)c2c(O)cc(O)cc2O)ccc1O',
                                     'name': 'homoeriodictyol chalcone',
                                     'reason': 'No phenyl acetate group found'},
                                 {   'smiles': 'O([C@H]1[C@H](O)[C@H](O[C@@H](O[C@H]2[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]2CO)O[C@@H]([C@H](O)[C@@H](NC(=O)C)CO)[C@H](O)CO)[C@H]1O)CO[C@H]3O[C@@H]([C@@H](O)[C@H](O)[C@@H]3O[C@@H]4O[C@@H]([C@@H](O[C@@H]5O[C@@H]([C@H](O)[C@H](O[C@]6(O[C@H]([C@H](NC(=O)C)[C@@H](O)C6)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]5O)CO)[C@H](O)[C@H]4NC(=O)C)CO)CO)[C@H]7O[C@@H]([C@@H](O)[C@H](O)[C@@H]7O[C@@H]8O[C@@H]([C@@H](O[C@@H]9O[C@@H]([C@H](O)[C@H](O)[C@H]9NC(=O)C)CO[C@]%10(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%10)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]8NC(=O)C)CO)CO',
                                     'name': '(2R,4S,5R,6R)-5-Acetamido-2-[[(2R,3R,4R,5R,6S)-5-acetamido-6-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2R,3S,4S,5S,6R)-2-[(2R,3R,4S,5S,6S)-2-[[(2S,3S,4S,5S,6R)-3-[(2S,3R,4R,5S,6R)-3-acetamido-5-[(2S,3R,4S,5S,6R)-4-[(2S,4S,5R,6R)-5-acetamido-2-carboxy-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxan-2-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-4-hydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]-6-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2R,3S,4R,5S)-5-acetamido-1,2,4,6-tetrahydroxyhexan-3-yl]oxy-4-hydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,5-dihydroxyoxan-4-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]oxy-4-hydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,4-dihydroxyoxan-2-yl]methoxy]-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'No phenyl acetate group found'},
                                 {   'smiles': 'O1C=2C(CC=C(C)C)=C(OC)C=C(O)C2C(=O)C=C1C3=CC=C(O)C=C3',
                                     'name': 'Artonin U',
                                     'reason': 'Contains excluded structural '
                                               'features'},
                                 {   'smiles': 'O=C(CCCCCCCCC)CCCCCC(O)=O',
                                     'name': '7-Keto palmitic acid',
                                     'reason': 'No phenyl acetate group found'},
                                 {   'smiles': 'CNC(=N)NCCC[C@H]([NH3+])C([O-])=O',
                                     'name': 'N(omega)-methyl-L-arginine '
                                             'zwitterion',
                                     'reason': 'No phenyl acetate group found'},
                                 {   'smiles': '[C@]([C@@](/C=C/CCCCCCCCCCC)(O)[H])(NC(=O)C(CCCCCCCCCCCCCCCCCCCCCC)O)([H])CO',
                                     'name': 'N-(2-hydroxy-tetracosanoyl)-hexadecasphing-4-enine',
                                     'reason': 'No phenyl acetate group found'},
                                 {   'smiles': 'COc1cc2oc(-c3ccc(O)cc3O)c(CC=C(C)C)c(=O)c2c(O)c1\\C=C\\C(C)C',
                                     'name': 'artocarpin',
                                     'reason': 'Contains excluded structural '
                                               'features'}],
    'sample_false_negatives': [   {   'smiles': 'CC(=O)Oc1ccc(cc1)[N+]([O-])=O',
                                      'name': '4-nitrophenyl acetate',
                                      'reason': 'Contains excluded structural '
                                                'features'}],
    'attempt': 3,
    'reasoning': "Let's analyze what went wrong and how to improve the "
                 'classifier:\n'
                 '\n'
                 '1. Main Issues from Results:\n'
                 '- High False Positive rate: The classifier incorrectly '
                 'identified many compounds as phenyl acetates\n'
                 '- Most false positives had acetate groups (-OC(=O)C) '
                 'attached to aromatic rings but were more complex structures\n'
                 "- The current pattern is too permissive and doesn't properly "
                 'exclude complex polycyclic compounds\n'
                 '\n'
                 '2. Key Improvements Needed:\n'
                 '\n'
                 'a) Pattern Refinement:\n'
                 '- Current pattern `c1([OH0]-[CX3](=[OX1])[CH3])ccccc1` is '
                 'too simple\n'
                 '- Need to better specify the phenol part (single benzene '
                 'ring)\n'
                 '- Need to better define the acetate connection\n'
                 '\n'
                 'b) Additional Checks:\n'
                 '- Need stricter verification of ring isolation\n'
                 '- Better handling of substituted phenyl rings\n'
                 '- More comprehensive negative patterns\n'
                 '\n'
                 '3. Proposed Solution:\n'
                 '\n'
                 "Here's an improved version:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 26,
    'num_false_positives': 71,
    'num_true_negatives': 142202,
    'num_false_negatives': 1,
    'num_negatives': None,
    'precision': 0.26804123711340205,
    'recall': 0.9629629629629629,
    'f1': 0.4193548387096774,
    'accuracy': 0.9994940267041462,
    'negative_predictive_value': 0.9999929677995542}