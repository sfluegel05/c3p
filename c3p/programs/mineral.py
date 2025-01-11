"""
Classifies: CHEBI:46662 mineral
"""
"""
Classifies: CHEBI:33686 mineral
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_mineral(smiles: str):
    """
    Determines if a molecule is a mineral based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a mineral, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get all atoms
    atoms = mol.GetAtoms()
    if len(atoms) == 0:
        return False, "No atoms found"
    
    # Count different types of atoms
    metal_ions = set()
    has_metal = False
    organic_atoms = 0
    inorganic_atoms = 0
    charged_atoms = 0
    
    # Define common metal atoms in minerals
    metals = {11: 'Na', 12: 'Mg', 13: 'Al', 19: 'K', 20: 'Ca', 24: 'Cr',
              25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn',
              38: 'Sr', 47: 'Ag', 48: 'Cd', 56: 'Ba', 57: 'La', 78: 'Pt',
              79: 'Au', 80: 'Hg', 82: 'Pb', 83: 'Bi', 55: 'Cs', 51: 'Sb',
              46: 'Pd', 33: 'As', 3: 'Li', 4: 'Be', 37: 'Rb', 81: 'Tl'}
    
    # Common anions/inorganic elements
    inorganic_elements = {8: 'O', 9: 'F', 15: 'P', 16: 'S', 17: 'Cl', 
                         35: 'Br', 53: 'I', 7: 'N', 5: 'B', 14: 'Si'}
    
    # Track heteroatoms and their charges
    heteroatom_charges = []
    aromatic_atoms = 0
    
    for atom in atoms:
        atomic_num = atom.GetAtomicNum()
        formal_charge = atom.GetFormalCharge()
        
        if atom.GetIsAromatic():
            aromatic_atoms += 1
            
        # Check for metal ions
        if atomic_num in metals:
            has_metal = True
            metal_ions.add(metals[atomic_num])
            if formal_charge != 0:
                charged_atoms += 1
                heteroatom_charges.append(formal_charge)
                
        # Count organic vs inorganic atoms
        elif atomic_num == 6:  # Carbon
            if len([n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]) > 1:
                organic_atoms += 1
        elif atomic_num in inorganic_elements:
            inorganic_atoms += 1
            if formal_charge != 0:
                charged_atoms += 1
                heteroatom_charges.append(formal_charge)
                
    # Check for water molecules
    water_pattern = Chem.MolFromSmiles("O")
    water_count = len(mol.GetSubstructMatches(water_pattern))
    
    # Check for complex organic groups
    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumRings()
    
    # Exclusion rules
    if aromatic_atoms > 0:
        return False, "Contains aromatic rings - likely organic"
    
    if ring_count > 2:
        return False, "Too many rings for a mineral"
        
    if organic_atoms > 4 and not any(x in ["Mg", "Ca", "Ba", "Sr"] for x in metal_ions):
        return False, "Too many C-C bonds for a mineral"
        
    # Must have either a metal or metalloid (As, Sb)
    if not has_metal and not any(x in metal_ions for x in ['As', 'Sb']):
        return False, "No metal or metalloid atoms found"
    
    # Simple metal halides and oxides
    simple_inorganic_pattern = all(a.GetAtomicNum() in list(metals.keys()) + [9, 17, 35, 53, 8] for a in atoms)
    if simple_inorganic_pattern:
        return True, f"Simple inorganic mineral containing {', '.join(metal_ions)}"
    
    # Check charge balance for ionic compounds
    if charged_atoms > 0:
        if sum(heteroatom_charges) != 0 and water_count == 0:
            return False, "Unbalanced charges in structure"
    
    # Classification logic
    if water_count > 0:
        if organic_atoms > 2 and not any(x in ["Mg", "Ca", "Ba", "Sr"] for x in metal_ions):
            return False, "Hydrated organic compound"
        return True, f"Hydrated mineral containing {', '.join(metal_ions)}"
    
    elif charged_atoms >= 2:
        if organic_atoms > 2 and not any(x in ["Mg", "Ca", "Ba", "Sr"] for x in metal_ions):
            return False, "Organic salt"
        return True, f"Ionic mineral containing {', '.join(metal_ions)}"
    
    elif simple_inorganic_pattern:
        return True, f"Covalent mineral containing {', '.join(metal_ions)}"
        
    return False, "Structure not consistent with mineral characteristics"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:46662',
                          'name': 'mineral',
                          'definition': 'In general, a mineral is a chemical '
                                        'substance that is normally '
                                        'crystalline formed and has been '
                                        'formed as a result of geological '
                                        'processes. The term also includes '
                                        'metamict substances (naturally '
                                        'occurring, formerly crystalline '
                                        'substances whose crystallinity has '
                                        'been destroyed by ionising radiation) '
                                        'and can include naturally occurring '
                                        'amorphous substances that have never '
                                        "been crystalline ('mineraloids') such "
                                        'as georgite and calciouranoite as '
                                        'well as substances formed by the '
                                        'action of geological processes on '
                                        "bigenic compounds ('biogenic "
                                        "minerals').",
                          'parents': ['CHEBI:59999'],
                          'xrefs': ['Wikipedia:Mineral'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.45502645502645495 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: [Ca++].[O-]S([O-])(=O)=O NAME: calcium '
               'sulfate REASON: CORRECT Hydrated mineral containing Ca\n'
               ' * SMILES: [Mg++].CCC([O-])=O.CCC([O-])=O NAME: magnesium '
               'dipropionate REASON: CORRECT Hydrated mineral containing Mg\n'
               ' * SMILES: [Ca++].[Ca++].[O-][Si]([O-])([O-])[O-] NAME: '
               'calcium silicate REASON: CORRECT Hydrated mineral containing '
               'Ca\n'
               ' * SMILES: [Ca++].[H]OP([O-])([O-])=O NAME: calcium '
               'hydrogenphosphate REASON: CORRECT Hydrated mineral containing '
               'Ca\n'
               ' * SMILES: [Ba++].[O-]C([O-])=O NAME: barium carbonate REASON: '
               'CORRECT Hydrated mineral containing Ba\n'
               ' * SMILES: O.O.Cl[Cu]Cl NAME: copper(II) chloride dihydrate '
               'REASON: CORRECT Hydrated mineral containing Cu\n'
               ' * SMILES: [Fe+3].[O-]P([O-])(=O)[O-] NAME: iron(3+) phosphate '
               'REASON: CORRECT Hydrated mineral containing Fe\n'
               ' * SMILES: [OH-].[OH-].[Ca++] NAME: calcium dihydroxide '
               'REASON: CORRECT Hydrated mineral containing Ca\n'
               ' * SMILES: [S--].[Fe+3].[As-] NAME: arsenopyrite REASON: '
               'CORRECT Ionic mineral containing As, Fe\n'
               ' * SMILES: [Cl-].[Cs+] NAME: caesium chloride REASON: CORRECT '
               'Ionic mineral containing Cs\n'
               ' * SMILES: [Fe++].[S-][S-] NAME: pyrite REASON: CORRECT Ionic '
               'mineral containing Fe\n'
               ' * SMILES: [Pd-2](Cl)(Cl)(Cl)(Cl)(Cl)Cl.[K+].[K+] NAME: '
               'Potassium hexachloropalladate(IV) REASON: CORRECT Ionic '
               'mineral containing Pd, K\n'
               ' * SMILES: [Ca+2].C(=O)([O-])[O-] NAME: calcium carbonate '
               'REASON: CORRECT Hydrated mineral containing Ca\n'
               ' * SMILES: [Cl-].[K+] NAME: potassium chloride REASON: CORRECT '
               'Ionic mineral containing K\n'
               ' * SMILES: Cl[Cu]Cl.O.O.O.O.O NAME: copper(II) chloride '
               'pentahydrate REASON: CORRECT Hydrated mineral containing Cu\n'
               ' * SMILES: [Ba++].CC([O-])=O.CC([O-])=O NAME: barium acetate '
               'REASON: CORRECT Hydrated mineral containing Ba\n'
               ' * SMILES: [Ba++].[O-][N+]([O-])=O.[O-][N+]([O-])=O NAME: '
               'barium nitrate REASON: CORRECT Hydrated mineral containing Ba\n'
               ' * SMILES: '
               'O.O.O.O.O.O.O.O.O.O.O.O.O.O.O.O.O.O.[Al+3].[Al+3].[O-]S([O-])(=O)=O.[O-]S([O-])(=O)=O.[O-]S([O-])(=O)=O '
               'NAME: aluminium sulfate octadecahydrate REASON: CORRECT '
               'Hydrated mineral containing Al\n'
               ' * SMILES: [S--].[S--].[Fe++].[Cu++] NAME: chalcopyrite '
               'REASON: CORRECT Ionic mineral containing Cu, Fe\n'
               ' * SMILES: [Mg++].[O-]S([O-])(=O)=O NAME: magnesium sulfate '
               'REASON: CORRECT Hydrated mineral containing Mg\n'
               ' * SMILES: P(=O)([O-])([O-])[O-].[Mg+2].[NH4+].O.O.O.O.O.O '
               'NAME: ammonium magnesium phosphate hexahydrate REASON: CORRECT '
               'Hydrated mineral containing Mg\n'
               ' * SMILES: [Cs+].[H]C([O-])=O NAME: caesium formate REASON: '
               'CORRECT Hydrated mineral containing Cs\n'
               ' * SMILES: '
               '[Na+].[Na+].[H]O[H].[H]O[H].[H]O[H].[H]O[H].[H]O[H].[H]O[H].[H]O[H].[H]O[H].[H]OB1O[B-]2(O[H])OB(O[H])O[B-](O[H])(O1)O2 '
               'NAME: disodium tetraborate octahydrate REASON: CORRECT '
               'Hydrated mineral containing Na\n'
               ' * SMILES: '
               '[Na+].[Na+].[H]OB1O[B-]2(O[H])OB(O[H])O[B-](O[H])(O1)O2 NAME: '
               'disodium tetraborate REASON: CORRECT Hydrated mineral '
               'containing Na\n'
               ' * SMILES: [O-]S([O-])(=O)=O.[Ba+2] NAME: barium sulfate '
               'REASON: CORRECT Hydrated mineral containing Ba\n'
               'False positives: SMILES: [Na+].[Na+].Cl[Pt--](Cl)(Cl)Cl NAME: '
               'sodium tetrachloroplatinate REASON: WRONGLY CLASSIFIED Ionic '
               'mineral containing Na, Pt\n'
               ' * SMILES: '
               '[Na+].[H][C@]12SC(C)(C)[C@@H](N1C(=O)[C@H]2NC(=O)c1c(OCC)ccc2ccccc12)C([O-])=O '
               'NAME: nafcillin sodium REASON: WRONGLY CLASSIFIED Metal '
               'organic salt containing Na\n'
               ' * SMILES: '
               '[K+].[H][C@]12SC(C)(C)[C@@H](N1C(=O)[C@H]2NC(=O)C(C)Oc1ccccc1)C([O-])=O '
               'NAME: phenethicillin potassium REASON: WRONGLY CLASSIFIED '
               'Metal organic salt containing K\n'
               ' * SMILES: [Na+].CC(O)C([O-])=O NAME: sodium lactate REASON: '
               'WRONGLY CLASSIFIED Hydrated mineral containing Na\n'
               ' * SMILES: '
               '[Na+].[H][C@]12SC(C)(C)[C@@H](N1C(=O)[C@H]2NC(=O)[C@]([H])(NC(=O)c1cnc2cccnc2c1O)c1ccccc1)C([O-])=O '
               'NAME: apalcillin sodium REASON: WRONGLY CLASSIFIED Metal '
               'organic salt containing Na\n'
               ' * SMILES: '
               'O.O.C/C=N/C(C(=O)[O-])C(C)O.[Fe+2].C/C=N/C(C(=O)[O-])C(C)O '
               'NAME: Ferrotrenine REASON: WRONGLY CLASSIFIED Metal organic '
               'salt containing Fe\n'
               ' * SMILES: '
               'P(=O)([O-])(OC[C@@H](O)C(=O)O[As](=O)([O-])[O-])[O-] NAME: '
               '1-arsono-3-phospho-D-glycerate(4-) REASON: WRONGLY CLASSIFIED '
               'Hydrated mineral containing As\n'
               ' * SMILES: O=N[Cr](N=O)(N=O)N=O NAME: tetranitrosylchromium '
               'REASON: WRONGLY CLASSIFIED Hydrated mineral containing Cr\n'
               ' * SMILES: '
               'CC1=C2N3[C@H]([C@H](CC([O-])=O)[C@@]2(C)CCC([O-])=O)[C@]2(C)[N+]4=C([C@@H](CCC([O-])=O)[C@]2(C)CC(N)=O)C(C)=C2[N+]5=C(C=C6[N+](=C1[C@@H](CCC([O-])=O)C6(C)C)[Co-3]345)[C@@H](CCC([O-])=O)[C@]2(C)CC(N)=O '
               'NAME: cob(I)yrinate a,c diamide(5-) REASON: WRONGLY CLASSIFIED '
               'Metal organic salt containing Co\n'
               ' * SMILES: [H][As]([O-])([O-])=O NAME: arsonate(2-) REASON: '
               'WRONGLY CLASSIFIED Hydrated mineral containing As\n'
               ' * SMILES: '
               '[Na+].[Na+].[H]O[H].[H]O[H].[H]O[H].[H]O[H].[H]O[H].[H]O[H].[H]O[H].[H]O[H].[H]O[H].[H]O[H].[O-]S([O-])(=O)=O '
               'NAME: sodium sulfate decahydrate REASON: WRONGLY CLASSIFIED '
               'Hydrated mineral containing Na\n'
               ' * SMILES: '
               'C1=2N3C(C=C4[N+]5=C(C=C6N7C(=CC8=[N+](C(=C1)C(=C8CCC([O-])=O)CO)[Fe-2]735)C(=C6C)CCC([O-])=O)C(=C4C)C=C)=C(C2C)[C@H](CC/C=C(\\C)/CC/C=C(/CCC=C(C)C)\\C)O '
               'NAME: Fe(II)-heme i(2-) REASON: WRONGLY CLASSIFIED Metal '
               'organic salt containing Fe\n'
               ' * SMILES: '
               '[Na+].[H]O[H].[H][C@]12SC(C)(C)[C@@H](N1C(=O)[C@H]2NC(=O)c1c(OC)cccc1OC)C([O-])=O '
               'NAME: methicillin sodium monohydrate REASON: WRONGLY '
               'CLASSIFIED Metal organic salt containing Na\n'
               ' * SMILES: [Na+].OCCN1CCN(CC1)CCS([O-])(=O)=O NAME: sodium '
               '2-[4-(2-hydroxyethyl)piperazin-1-yl]ethanesulfonate REASON: '
               'WRONGLY CLASSIFIED Hydrated mineral containing Na\n'
               ' * SMILES: Br[Ag--](Br)Br NAME: tribromoargentate(2-) REASON: '
               'WRONGLY CLASSIFIED Mineral containing Ag\n'
               ' * SMILES: [S]12[Fe]3[S]4[Fe]1[S]1[Fe+]2[S]3[Fe+]41 NAME: '
               'tetra-mu3-sulfido-tetrairon(2+) REASON: WRONGLY CLASSIFIED '
               'Ionic mineral containing Fe\n'
               ' * SMILES: [Cs+].[Cs+].Br[Cr--](Br)(Br)Br NAME: dicaesium '
               'tetrabromochromate REASON: WRONGLY CLASSIFIED Ionic mineral '
               'containing Cs, Cr\n'
               ' * SMILES: [O-][Mn](=O)(=O)=O NAME: permanganate REASON: '
               'WRONGLY CLASSIFIED Hydrated mineral containing Mn\n'
               ' * SMILES: [H][N]([H])([H])[Hg++][N]([H])([H])[H] NAME: '
               'diamminemercury(2+) REASON: WRONGLY CLASSIFIED Mineral '
               'containing Hg\n'
               ' * SMILES: [K+].OCCN(CCO)C([S-])=S NAME: Potassium '
               'bis(2-hydroxyethyl)dithiocarbamate REASON: WRONGLY CLASSIFIED '
               'Hydrated mineral containing K\n'
               ' * SMILES: N#C[Ni-4](C#N)(C#N)C#N NAME: '
               'tetracyanonickelate(4-) REASON: WRONGLY CLASSIFIED Mineral '
               'containing Ni\n'
               ' * SMILES: I[Hg--](I)(I)I NAME: tetraiodomercurate(2-) REASON: '
               'WRONGLY CLASSIFIED Mineral containing Hg\n'
               ' * SMILES: '
               '[NH4+].[NH4+].[NH4+].N#C[Fe-3](C#N)(C#N)(C#N)(C#N)C#N NAME: '
               'ammonium hexacyanoferrate(3-) REASON: WRONGLY CLASSIFIED Ionic '
               'mineral containing Fe\n'
               ' * SMILES: F[Cr-](F)(F)(F)(F)F NAME: hexafluorochromate(1-) '
               'REASON: WRONGLY CLASSIFIED Mineral containing Cr\n'
               ' * SMILES: [NH2-].[Na+] NAME: sodium amide REASON: WRONGLY '
               'CLASSIFIED Ionic mineral containing Na\n'
               'False negatives: SMILES: Cl[La](Cl)Cl NAME: lanthanum '
               'trichloride REASON: MISSED No ionic character or coordination '
               'bonds found\n'
               ' * SMILES: [Sb](F)(F)(F)(F)F NAME: antimony pentafluoride '
               'REASON: MISSED No ionic character or coordination bonds found\n'
               ' * SMILES: [Ni]=S=[Ni]=S=[Ni] NAME: heazlewoodite REASON: '
               'MISSED No ionic character or coordination bonds found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C1C=C(N=C2N1NN=N2)C3=CC=CC=C3',
                                     'name': '5-phenyl-1,7-dihydrotetrazolo[1,5-a]pyrimidine',
                                     'reason': 'Contains aromatic rings - '
                                               'likely organic'},
                                 {   'smiles': 'O(C(=O)CCCCCCCCCCC/C=C\\CCCCCCCC)[C@H](COC(=O)CCCCCCCCC/C=C\\CCCCCC)COC(=O)CCCCCCC/C=C\\CCCC',
                                     'name': 'TG(14:1(9Z)/22:1(13Z)/18:1(11Z))',
                                     'reason': 'Too many C-C bonds for a '
                                               'mineral'},
                                 {   'smiles': 'CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@@H](O[C@@H]3O[C@H](CO[C@H]4O[C@H](CO[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H](O)[C@H](O[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H]4O)[C@@H](O)[C@H](O[C@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]4O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4NC(C)=O)[C@@H]3O)[C@H](O)[C@H]2NC(C)=O)[C@@H]1O',
                                     'name': 'beta-D-GlcpNAc-(1->2)-alpha-D-Manp-(1->3)-{alpha-D-Manp-(1->3)-[alpha-D-Manp-(1->6)]-alpha-D-Manp-(1->6)}-beta-D-Manp-(1->4)-beta-GlcpNAc-(1->4)-beta-D-GlcpNAc',
                                     'reason': 'Too many rings for a mineral'},
                                 {   'smiles': '[C@@]12([C@@]([C@]([C@@H](CC1)C)(CCC(CC)C)C)(CCC[C@@H]2C)[H])C',
                                     'name': 'clerodane',
                                     'reason': 'Too many C-C bonds for a '
                                               'mineral'},
                                 {   'smiles': 'S(O[C@@H]1[C@H](O)[C@@H](O)[C@H](O[C@@H]([C@@H](O)[C@H](O)CO[C@]2(O[C@H]([C@H](NC(=O)C)[C@@H](O)C2)[C@H](O)[C@H](O)CO)C(O)=O)[C@@H](NC(=O)C)CO)O[C@@H]1CO)(O)(=O)=O',
                                     'name': '(2R,4S,5R,6R)-5-Acetamido-2-[(2R,3S,4R,5S)-5-acetamido-4-[(2R,3R,4R,5R,6R)-3,4-dihydroxy-6-(hydroxymethyl)-5-sulfooxyoxan-2-yl]oxy-2,3,6-trihydroxyhexoxy]-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'Too many C-C bonds for a '
                                               'mineral'},
                                 {   'smiles': 'O=C(NC(CC(C)C)C(=O)N[C@@H](CCCN=C(N)N)C=O)C(NC(=O)CC)CC(C)C',
                                     'name': 'Leupeptin Pr-LL',
                                     'reason': 'Too many C-C bonds for a '
                                               'mineral'},
                                 {   'smiles': 'COC1=CC=C(C=C1)CCNC(=O)C(=C2C3=CC=CC=C3C(=N2)NC(=O)C4=CC=CS4)C#N',
                                     'name': 'N-[3-[1-cyano-2-[2-(4-methoxyphenyl)ethylamino]-2-oxoethylidene]-1-isoindolyl]-2-thiophenecarboxamide',
                                     'reason': 'Contains aromatic rings - '
                                               'likely organic'},
                                 {   'smiles': 'O=C1O[C@@H](CC=2C1=C(OC)C(O)=CC2)CCC[C@@H](O)C',
                                     'name': 'Penicimarin C',
                                     'reason': 'Contains aromatic rings - '
                                               'likely organic'},
                                 {   'smiles': 'CO[C@@H]1[C@@H]2[C@@H](C[C@@H]3[C@@](O2)([C@H]([C@@H]([C@H](O3)C(=O)OC)O)O)O)OC1N4C=NC5=C4N=CN=C5N',
                                     'name': 'LSM-4497',
                                     'reason': 'Contains aromatic rings - '
                                               'likely organic'},
                                 {   'smiles': 'C[C@@H]1[C@H](O)CCC2=CC[C@H](C[C@]12C)C(C)=C',
                                     'name': '1-deoxycapsidiol',
                                     'reason': 'Too many C-C bonds for a '
                                               'mineral'}],
    'sample_false_negatives': [   {   'smiles': '[Ni]=S=[Ni]=S=[Ni]',
                                      'name': 'heazlewoodite',
                                      'reason': 'Structure not consistent with '
                                                'mineral characteristics'}],
    'attempt': 1,
    'reasoning': 'Let me analyze the previous outcomes and suggest '
                 'improvements:\n'
                 '\n'
                 '1. False Positives Analysis:\n'
                 '- Many organic metal salts (like nafcillin sodium) were '
                 'incorrectly classified as minerals\n'
                 '- Metal-organic complexes (like Fe-heme) were wrongly '
                 'classified\n'
                 '- Simple metal coordination compounds were incorrectly '
                 'included\n'
                 '- The program was too lenient with organic content\n'
                 '\n'
                 '2. False Negatives Analysis:\n'
                 '- Simple metal halides (like LaCl3) were missed due to lack '
                 'of explicit charges\n'
                 '- Binary compounds (like SbF5) were missed due to lack of '
                 'ionic character\n'
                 '- Some mineral structures with covalent bonds (like '
                 'heazlewoodite) were missed\n'
                 '\n'
                 '3. Key Issues to Address:\n'
                 '- Need stricter rules for organic content\n'
                 '- Need to handle both ionic and covalent mineral structures\n'
                 '- Need better handling of coordination compounds\n'
                 '- Need to exclude biological metal complexes\n'
                 '- Need to better detect simple metal halides and oxides\n'
                 '\n'
                 "Here's my improved approach:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 45,
    'num_false_positives': 100,
    'num_true_negatives': 20508,
    'num_false_negatives': 1,
    'num_negatives': None,
    'precision': 0.3103448275862069,
    'recall': 0.9782608695652174,
    'f1': 0.4712041884816754,
    'accuracy': 0.9951099060714631,
    'negative_predictive_value': 0.9999512409186211}