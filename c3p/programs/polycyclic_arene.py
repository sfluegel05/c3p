"""
Classifies: CHEBI:33848 polycyclic arene
"""
"""
Classifies: CHEBI:33853 polycyclic arene
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polycyclic_arene(smiles: str):
    """
    Determines if a molecule is a polycyclic arene based on its SMILES string.
    A polycyclic arene is a polycyclic aromatic hydrocarbon consisting only of carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polycyclic arene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if all atoms are carbon (polycyclic arenes are pure hydrocarbons)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            return False, "Contains non-carbon atoms"

    # Get all rings
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # Check if there are at least 2 rings
    if len(rings) < 2:
        return False, "Not enough rings to be polycyclic"

    # Check if rings are fused and aromatic
    aromatic_rings = 0
    for ring in rings:
        # Check if all atoms in the ring are aromatic
        if all(mol.GetAtomWithIdx(atom).GetIsAromatic() for atom in ring):
            aromatic_rings += 1
            # Check if this ring shares at least 2 atoms with another ring
            shared_atoms = 0
            for other_ring in rings:
                if ring != other_ring:
                    shared_atoms = len(set(ring).intersection(other_ring))
                    if shared_atoms >= 2:
                        break
            else:
                return False, "Rings are not properly fused"

    # Need at least 2 aromatic rings
    if aromatic_rings < 2:
        return False, f"Only {aromatic_rings} aromatic ring(s) found, need at least 2"

    return True, "Contains multiple fused aromatic rings consisting only of carbon atoms"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33848',
                          'name': 'polycyclic arene',
                          'definition': 'A polycyclic aromatic hydrocarbon.',
                          'parents': ['CHEBI:33658', 'CHEBI:33666'],
                          'xrefs': [   'PMID:15198916',
                                       'PMID:25679824',
                                       'Wikipedia:Polycyclic_aromatic_hydrocarbon'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.6062992125984252 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: Oc1ccc2cc3ccc4cccc5ccc(c2c1O)c3c45 '
               'NAME: 9,10-Dihydroxybenzo[a]pyrene REASON: CORRECT Contains '
               'multiple fused aromatic rings consisting only of carbon atoms\n'
               ' * SMILES: [O-][N+](=O)c1c2ccccc2c2ccc3cccc4ccc1c2c34 NAME: '
               '6-Nitrobenzo[a]pyrene REASON: CORRECT Contains multiple fused '
               'aromatic rings consisting only of carbon atoms\n'
               ' * SMILES: '
               'c1cc2cccc3c4ccc5c6ccc7c8cccc9cccc(c%10ccc(c%11ccc(c(c1)c23)c4c5%11)c6c7%10)c89 '
               'NAME: quaterrylene REASON: CORRECT Contains multiple fused '
               'aromatic rings consisting only of carbon atoms\n'
               ' * SMILES: c1ccc2ccccc2c1 NAME: naphthalene REASON: CORRECT '
               'Contains multiple fused aromatic rings consisting only of '
               'carbon atoms\n'
               ' * SMILES: Cc1c2ccccc2c(C)c2c1ccc1ccccc21 NAME: '
               '7,12-dimethyltetraphene REASON: CORRECT Contains multiple '
               'fused aromatic rings consisting only of carbon atoms\n'
               ' * SMILES: '
               'c1cc2ccc3c4ccc5ccc6ccc7ccc8c9ccc%10ccc1c1c2c3c(c9c%101)c1c4c5c6c7c81 '
               'NAME: circumbiphenyl REASON: CORRECT Contains multiple fused '
               'aromatic rings consisting only of carbon atoms\n'
               ' * SMILES: C1C=C2C(=C1)c1cccc3cccc2c13 NAME: '
               '8H-cyclopenta[a]acenaphthylene REASON: CORRECT Contains '
               'multiple fused aromatic rings consisting only of carbon atoms\n'
               ' * SMILES: Oc1cc2ccc3cc4ccccc4c4ccc(c1)c2c34 NAME: '
               '2-Hydroxybenzo[a]pyrene REASON: CORRECT Contains multiple '
               'fused aromatic rings consisting only of carbon atoms\n'
               ' * SMILES: c1ccc2c(c1)cc1ccc3cccc4c5ccccc5c2c1c34 NAME: '
               'dibenzo[a,l]pyrene REASON: CORRECT Contains multiple fused '
               'aromatic rings consisting only of carbon atoms\n'
               ' * SMILES: c1cc2cc3cccc4c5cccc6cc7cccc8c(c1)c2c(c34)c(c78)c56 '
               'NAME: bisanthene REASON: CORRECT Contains multiple fused '
               'aromatic rings consisting only of carbon atoms\n'
               ' * SMILES: Oc1c2ccccc2c2ccc3cccc4ccc1c2c34 NAME: '
               '6-Hydroxybenzo[a]pyrene REASON: CORRECT Contains multiple '
               'fused aromatic rings consisting only of carbon atoms\n'
               ' * SMILES: c1cc2cccc3c4cccc5cccc(c(c1)c23)c45 NAME: perylene '
               'REASON: CORRECT Contains multiple fused aromatic rings '
               'consisting only of carbon atoms\n'
               ' * SMILES: c1ccc2c(c1)cc1c3ccccc3c3cc4ccccc4c4ccc2c1c34 NAME: '
               'Dibenzo[h,rst]pentaphene REASON: CORRECT Contains multiple '
               'fused aromatic rings consisting only of carbon atoms\n'
               ' * SMILES: c1ccc2c(c1)c1cccc3c4ccccc4c4cccc2c4c13 NAME: '
               'Dibenzo[e,l]pyrene REASON: CORRECT Contains multiple fused '
               'aromatic rings consisting only of carbon atoms\n'
               ' * SMILES: Oc1c(O)c2c3ccccc3cc3ccc4cccc1c4c23 NAME: '
               '11,12-Dihydroxybenzo[a]pyrene REASON: CORRECT Contains '
               'multiple fused aromatic rings consisting only of carbon atoms\n'
               ' * SMILES: c1ccc2c(c1)cc1ccc3cccc4ccc2c1c34 NAME: '
               'benzo[a]pyrene REASON: CORRECT Contains multiple fused '
               'aromatic rings consisting only of carbon atoms\n'
               ' * SMILES: c1ccc-2c(c1)-c1ccccc1-c1ccccc1-c1ccccc-21 NAME: '
               'tetraphenylene REASON: CORRECT Contains multiple fused '
               'aromatic rings consisting only of carbon atoms\n'
               ' * SMILES: '
               'c1ccc2c(c1)ccc1ccc3ccc4ccc5ccc6ccc7ccc8ccccc8c7c6c5c4c3c21 '
               'NAME: nonahelicene REASON: CORRECT Contains multiple fused '
               'aromatic rings consisting only of carbon atoms\n'
               ' * SMILES: c1ccc2cc3c4cccc5ccc6cccc(c3cc2c1)c6c45 NAME: '
               'Naphtho[2,3-e]pyrene REASON: CORRECT Contains multiple fused '
               'aromatic rings consisting only of carbon atoms\n'
               ' * SMILES: c1ccc2cc3cccc4cccc(cc2c1)c34 NAME: pleiadene '
               'REASON: CORRECT Contains multiple fused aromatic rings '
               'consisting only of carbon atoms\n'
               ' * SMILES: '
               'c1ccc2c(c1)ccc1ccc3ccc4ccc5ccc6ccc7ccccc7c6c5c4c3c21 NAME: '
               'octahelicene REASON: CORRECT Contains multiple fused aromatic '
               'rings consisting only of carbon atoms\n'
               ' * SMILES: c1cc2ccc3cc4cccc5ccc6cc(c1)c2c3c6c45 NAME: '
               'Anthanthrene REASON: CORRECT Contains multiple fused aromatic '
               'rings consisting only of carbon atoms\n'
               ' * SMILES: '
               'c1cc2cc3cc4ccc5cc6cc7ccc8cc9ccc%10cc%11cc%12ccc%13cc%14cc%15ccc%16cc1c1c2c2c3c3c4c5c4c6c5c7c8c6c9c%10c7c%11c8c%12c%13c9c%14c%10c%15c%16c1c1c2c2c3c4c3c5c6c7c4c8c9c(c%101)c2c34 '
               'NAME: circumcircumpyrene REASON: CORRECT Contains multiple '
               'fused aromatic rings consisting only of carbon atoms\n'
               ' * SMILES: '
               'c1cc2c3ccc4c5ccc6c7ccc8c9ccc%10c%11ccc%12c1c2c1c(c%11%12)c(c9%10)c(c78)c(c56)c1c34 '
               'NAME: hexacyclopentacoronene REASON: CORRECT Contains multiple '
               'fused aromatic rings consisting only of carbon atoms\n'
               ' * SMILES: c1ccc2c(c1)c1cccc3c4c5ccccc5c5cccc(c2c13)c45 NAME: '
               'rubicene REASON: CORRECT Contains multiple fused aromatic '
               'rings consisting only of carbon atoms\n'
               'False positives: SMILES: '
               'O1C(CCC=2C1=CC(OC3OC(C(O)C(O)C3O)C(O)=O)=C(C2OC)C=4C(=O)C=5C(OC4)=CC(O)=C(O)C5)(C)C '
               'NAME: '
               '6-{[6-(6,7-dihydroxy-4-oxo-4H-chromen-3-yl)-5-methoxy-2,2-dimethyl-3,4-dihydro-2H-1-benzopyran-7-yl]oxy}-3,4,5-trihydroxyoxane-2-carboxylic '
               'acid REASON: WRONGLY CLASSIFIED Contains multiple fused '
               'aromatic rings consisting only of carbon atoms\n'
               ' * SMILES: '
               'O=C(O)[C@]1([C@H]2[C@@](OC=3C=C4C5=C(C=CC=C5)NC4=CC3CC2)(CC[C@@H]1O)C)C '
               'NAME: Oxiamycin REASON: WRONGLY CLASSIFIED Contains multiple '
               'fused aromatic rings consisting only of carbon atoms\n'
               ' * SMILES: '
               'C1CCC(C1)CC#CC2=CC=C(C=C2)[C@@H]3[C@@H]4CN(CC(=O)N4[C@@H]3CO)C(=O)C5=CC=C(C=C5)F '
               'NAME: '
               '(6R,7R,8S)-7-[4-(3-cyclopentylprop-1-ynyl)phenyl]-4-[(4-fluorophenyl)-oxomethyl]-8-(hydroxymethyl)-1,4-diazabicyclo[4.2.0]octan-2-one '
               'REASON: WRONGLY CLASSIFIED Contains multiple fused aromatic '
               'rings consisting only of carbon atoms\n'
               ' * SMILES: '
               'CC=CC1=CC=C(C=C1)[C@@H]2[C@@H](N(C23CN(C3)C(=O)C)C(=O)C4=NC=CN=C4)CO '
               'NAME: LSM-39875 REASON: WRONGLY CLASSIFIED Contains multiple '
               'fused aromatic rings consisting only of carbon atoms\n'
               ' * SMILES: CN1C2=CC=CC=C2OC1=S NAME: '
               '3-methyl-1,3-benzoxazole-2-thione REASON: WRONGLY CLASSIFIED '
               'Contains multiple fused aromatic rings consisting only of '
               'carbon atoms\n'
               ' * SMILES: '
               'COC1=CC=CC(=C1)NC(=O)N2C[C@@H]3[C@H]([C@@H](N3C(=O)C2)CO)C4=CC=C(C=C4)C#CCC5CCCC5 '
               'NAME: '
               '(6S,7R,8R)-7-[4-(3-cyclopentylprop-1-ynyl)phenyl]-8-(hydroxymethyl)-N-(3-methoxyphenyl)-2-oxo-1,4-diazabicyclo[4.2.0]octane-4-carboxamide '
               'REASON: WRONGLY CLASSIFIED Contains multiple fused aromatic '
               'rings consisting only of carbon atoms\n'
               ' * SMILES: '
               'C1=CC=C(C=C1)C2=C(N3C(=O)C(=CC4=CC(=CC=C4)O)SC3=N2)C5=CC=CC=C5 '
               'NAME: '
               '2-[(3-hydroxyphenyl)methylidene]-5,6-diphenyl-3-imidazo[2,1-b]thiazolone '
               'REASON: WRONGLY CLASSIFIED Contains multiple fused aromatic '
               'rings consisting only of carbon atoms\n'
               ' * SMILES: '
               'O=C1C(N=C2C(=O)C=3C=4C(=C(C)C=C(C4C2=O)O)C=5O[C@@H](C(C5C3O)(C)C)C)=C(O)C=6C(O)=CC(=C7C6C1=C(O)C=8C([C@@H](C)OC78)(C)C)C '
               'NAME: Phenalenone derivative REASON: WRONGLY CLASSIFIED '
               'Contains multiple fused aromatic rings consisting only of '
               'carbon atoms\n'
               ' * SMILES: '
               'ClC1=C(OC)C=CC(=C1)C[C@H](NC(=O)/C=C/C[C@H](O)[C@@H]([C@@H](O)[C@@H](O)C2=CC=CC=C2)C)C(=O)O '
               'NAME: Cryptophycin G REASON: WRONGLY CLASSIFIED Contains '
               'multiple fused aromatic rings consisting only of carbon atoms\n'
               ' * SMILES: '
               'CC1=CC(=C(C=C1)OC)NC(=O)C2=CC3=CC=CC=C3C=C2OCC(=O)NCC4=CC=CO4 '
               'NAME: '
               '3-[2-(2-furanylmethylamino)-2-oxoethoxy]-N-(2-methoxy-5-methylphenyl)-2-naphthalenecarboxamide '
               'REASON: WRONGLY CLASSIFIED Contains multiple fused aromatic '
               'rings consisting only of carbon atoms\n'
               ' * SMILES: C1=2C(O[C@@H](CC1=O)C3=CC=C(C=C3)OC)=CC(=C(C2O)O)OC '
               "NAME: (2S)-carthamidin-4',7-dimethyl ether REASON: WRONGLY "
               'CLASSIFIED Contains multiple fused aromatic rings consisting '
               'only of carbon atoms\n'
               ' * SMILES: '
               'O1C(C(O)C(O)C(O)C1OC2=CC3=C(C=C2)C=C(OC3=O)/C=C/CC)C(O)=O '
               'NAME: '
               '6-({3-[(1E)-but-1-en-1-yl]-1-oxo-1H-isochromen-7-yl}oxy)-3,4,5-trihydroxyoxane-2-carboxylic '
               'acid REASON: WRONGLY CLASSIFIED Contains multiple fused '
               'aromatic rings consisting only of carbon atoms\n'
               ' * SMILES: '
               'O(C(C1OC2=C(C1)C(O)=C3C(OC(=O)C=C3)=C2)(C)C)C4OC(C(O)C(O)C4O)CO '
               'NAME: Celereoside REASON: WRONGLY CLASSIFIED Contains multiple '
               'fused aromatic rings consisting only of carbon atoms\n'
               ' * SMILES: '
               '[Cl-].CN(C)c1ccc(cc1)C(=C1C=CC(C=C1)=[N+](C)C)c1ccc(Nc2ccccc2)c2ccccc12 '
               'NAME: victoria blue B REASON: WRONGLY CLASSIFIED Contains '
               'multiple fused aromatic rings consisting only of carbon atoms\n'
               ' * SMILES: '
               'C[C@H]1CN(C(=O)CC2=C(C=CC(=C2)NC(=O)CN3C=NN=N3)O[C@H]1CN(C)CC4=CC=C(C=C4)OC5=CC=CC=C5)[C@@H](C)CO '
               'NAME: '
               'N-[(2R,3S)-5-[(2S)-1-hydroxypropan-2-yl]-3-methyl-2-[[methyl-[(4-phenoxyphenyl)methyl]amino]methyl]-6-oxo-2,3,4,7-tetrahydro-1,5-benzoxazonin-9-yl]-2-(1-tetrazolyl)acetamide '
               'REASON: WRONGLY CLASSIFIED Contains multiple fused aromatic '
               'rings consisting only of carbon atoms\n'
               ' * SMILES: '
               'CC[C@@]1(O)C[C@H](O[C@H]2C[C@@H]([C@H](O)[C@H](C)O2)N(C)C)c2c(O)c3C(=O)c4c(O)cccc4C(=O)c3c(O)c2[C@H]1C(=O)OC '
               'NAME: epsilon-Rhodomycin T REASON: WRONGLY CLASSIFIED Contains '
               'multiple fused aromatic rings consisting only of carbon atoms\n'
               ' * SMILES: '
               'CN1C(=NN=C1SCC(=O)N2C3=C(CCCC3)C4=CC=CC=C42)C5=CC=CO5 NAME: '
               '2-[[5-(2-furanyl)-4-methyl-1,2,4-triazol-3-yl]thio]-1-(1,2,3,4-tetrahydrocarbazol-9-yl)ethanone '
               'REASON: WRONGLY CLASSIFIED Contains multiple fused aromatic '
               'rings consisting only of carbon atoms\n'
               ' * SMILES: COc1ccc(CCCN2CCC(COC(c3ccccc3)c3ccccc3)CC2)cc1 '
               'NAME: UK 78282 REASON: WRONGLY CLASSIFIED Contains multiple '
               'fused aromatic rings consisting only of carbon atoms\n'
               ' * SMILES: '
               'C1=2[C@@H](C(=C(OC1=CC(=CC2O)O)C3=CC=C(C=C3)O)[O-])O NAME: '
               '(4S)-2,3-dehydroleucopelargonidin(1-) REASON: WRONGLY '
               'CLASSIFIED Contains multiple fused aromatic rings consisting '
               'only of carbon atoms\n'
               ' * SMILES: '
               'O=C1N[C@H](C(=O)N[C@H](C(=O)N([C@@H](C(NC2=C1C=CC=C2)=O)[C@@H](O)C3=CC=CC=C3)C)C(C)C)C '
               'NAME: Asperterrestride B REASON: WRONGLY CLASSIFIED Contains '
               'multiple fused aromatic rings consisting only of carbon atoms\n'
               ' * SMILES: '
               'C[C@@H]1CCCCO[C@H]([C@@H](CN(C(=O)C2=C(O1)C=CC(=C2)NC(=O)C3=CC=NC=C3)[C@H](C)CO)C)CN(C)CC4=CC=NC=C4 '
               'NAME: '
               'N-[(3R,9R,10R)-12-[(2R)-1-hydroxypropan-2-yl]-3,10-dimethyl-9-[[methyl(pyridin-4-ylmethyl)amino]methyl]-13-oxo-2,8-dioxa-12-azabicyclo[12.4.0]octadeca-1(14),15,17-trien-16-yl]-4-pyridinecarboxamide '
               'REASON: WRONGLY CLASSIFIED Contains multiple fused aromatic '
               'rings consisting only of carbon atoms\n'
               ' * SMILES: '
               'O(C=1C=C(C=CC1O)/C=C/C(O)=CC(=O)\\C=C\\C2=CC=CC=C2)C NAME: '
               '1-(4-Hydroxy-3-methoxyphenyl)-7-phenyl-5-hydroxy-1,4,6-heptatriene-3-one '
               'REASON: WRONGLY CLASSIFIED Contains multiple fused aromatic '
               'rings consisting only of carbon atoms\n'
               ' * SMILES: '
               'CC1=CC(=CC(=C1)OCC(=O)OCC(=O)C(C#N)C2=NC3=CC=CC=C3N2)C NAME: '
               '2-(3,5-dimethylphenoxy)acetic acid '
               '[3-(1H-benzimidazol-2-yl)-3-cyano-2-oxopropyl] ester REASON: '
               'WRONGLY CLASSIFIED Contains multiple fused aromatic rings '
               'consisting only of carbon atoms\n'
               ' * SMILES: '
               'Cl[C@H]1[C@](C=C)([C@@H]2NC(=O)S[C@@]23C=4C=5C(=CC=CC5C([C@H]3C1)(C)C)NC4)C '
               'NAME: Hapalindole T REASON: WRONGLY CLASSIFIED Contains '
               'multiple fused aromatic rings consisting only of carbon atoms\n'
               ' * SMILES: '
               'CCS(=O)(=O)N1CC2(C1)[C@@H]([C@@H](N2C(=O)C3=CC=CC=C3)CO)C4=CC=CC=C4 '
               'NAME: '
               '[(2R,3R)-6-ethylsulfonyl-2-(hydroxymethyl)-3-phenyl-1,6-diazaspiro[3.3]heptan-1-yl]-phenylmethanone '
               'REASON: WRONGLY CLASSIFIED Contains multiple fused aromatic '
               'rings consisting only of carbon atoms\n'
               'False negatives: NONE\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'OC[C@H]1O[C@H](C[C@@H]1O)N1C=NC2=C1N=CNC[C@H]2O',
                                     'name': 'pentostatin',
                                     'reason': 'Contains non-carbon atoms'},
                                 {   'smiles': 'O1C(CCC=2C1=CC(OC3OC(C(O)C(O)C3O)C(O)=O)=C(C2OC)C=4C(=O)C=5C(OC4)=CC(O)=C(O)C5)(C)C',
                                     'name': '6-{[6-(6,7-dihydroxy-4-oxo-4H-chromen-3-yl)-5-methoxy-2,2-dimethyl-3,4-dihydro-2H-1-benzopyran-7-yl]oxy}-3,4,5-trihydroxyoxane-2-carboxylic '
                                             'acid',
                                     'reason': 'Contains non-carbon atoms'},
                                 {   'smiles': 'O=C1O[C@@H]([C@@H](OC)C=CC=C(C[C@@H](C)[C@@H]([C@@H]([C@@H]([C@@H](C=C(C=C1OC)C)C)O)C)O)C)[C@H]([C@@H](O)[C@@H]([C@@]2(O[C@H](/C=C/C)[C@@H](C)[C@@H](C2)OC3OC(C(O)C(C3)O)C)O)C)C',
                                     'name': 'Concanamycin D',
                                     'reason': 'Contains non-carbon atoms'},
                                 {   'smiles': 'O1[C@@H](O[C@@H]2[C@@H](O)[C@@H](O[C@@H]([C@@H]2O)CO)O[C@H]3[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]3CO)O)[C@H](NC(=O)C)[C@@H](O[C@@H]4O[C@@H]([C@H](O)[C@H](O)[C@H]4O)CO)[C@H](O[C@@H]5O[C@H]([C@@H](O)[C@@H](O)[C@@H]5O)C)[C@H]1CO',
                                     'name': 'N-[(3R,4R,5S,6R)-5-[(2S,3R,4S,5S,6R)-4-[(2S,3R,4R,5S,6R)-3-Acetamido-6-(hydroxymethyl)-4-[(2R,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-5-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-2-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-2,4-dihydroxy-6-(hydroxymethyl)oxan-3-yl]acetamide',
                                     'reason': 'Contains non-carbon atoms'},
                                 {   'smiles': 'O=C(O)[C@]1([C@H]2[C@@](OC=3C=C4C5=C(C=CC=C5)NC4=CC3CC2)(CC[C@@H]1O)C)C',
                                     'name': 'Oxiamycin',
                                     'reason': 'Contains non-carbon atoms'},
                                 {   'smiles': 'C1CCC(C1)CC#CC2=CC=C(C=C2)[C@@H]3[C@@H]4CN(CC(=O)N4[C@@H]3CO)C(=O)C5=CC=C(C=C5)F',
                                     'name': '(6R,7R,8S)-7-[4-(3-cyclopentylprop-1-ynyl)phenyl]-4-[(4-fluorophenyl)-oxomethyl]-8-(hydroxymethyl)-1,4-diazabicyclo[4.2.0]octan-2-one',
                                     'reason': 'Contains non-carbon atoms'},
                                 {   'smiles': 'S(=O)(C(SSCCC)CC)CCC',
                                     'name': 'Propyl 1-(propylsulfinyl)propyl '
                                             'disulfide',
                                     'reason': 'Contains non-carbon atoms'},
                                 {   'smiles': 'ClC=1C(=O)[C@@H]([C@@](O)(C/C=C\\CCCCC)C1)C[C@H](OC(=O)C)[C@@H](OC(=O)C)CCCC(OC)=O',
                                     'name': 'punaglandin 6',
                                     'reason': 'Contains non-carbon atoms'},
                                 {   'smiles': 'CCCCCCCCCCCCCCCCCCCCCCCC[C@H](O)C(=O)N[C@@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)[C@H](O)[C@@H](O)CCCCCCCCCCCCCC',
                                     'name': 'Ins-1-P-Cer(t18:0/2-OH-26:0)',
                                     'reason': 'Contains non-carbon atoms'},
                                 {   'smiles': 'CCC(C)(C)C(=O)OC1CC(C=C2C1[C@H]([C@H](C=C2)C)CC[C@@H]3CC(CC(=O)O3)O)C',
                                     'name': '2,2-dimethylbutanoic acid '
                                             '[(7S,8S)-8-[2-[(2R)-4-hydroxy-6-oxo-2-oxanyl]ethyl]-3,7-dimethyl-1,2,3,7,8,8a-hexahydronaphthalen-1-yl] '
                                             'ester',
                                     'reason': 'Contains non-carbon atoms'}],
    'sample_false_negatives': [   {   'smiles': 'Oc1ccc2cc3ccc4cccc5ccc(c2c1O)c3c45',
                                      'name': '9,10-Dihydroxybenzo[a]pyrene',
                                      'reason': 'Contains non-carbon atoms'},
                                  {   'smiles': '[O-][N+](=O)c1c2ccccc2c2ccc3cccc4ccc1c2c34',
                                      'name': '6-Nitrobenzo[a]pyrene',
                                      'reason': 'Contains non-carbon atoms'},
                                  {   'smiles': 'Oc1cc2ccc3cc4ccccc4c4ccc(c1)c2c34',
                                      'name': '2-Hydroxybenzo[a]pyrene',
                                      'reason': 'Contains non-carbon atoms'},
                                  {   'smiles': 'Oc1c2ccccc2c2ccc3cccc4ccc1c2c34',
                                      'name': '6-Hydroxybenzo[a]pyrene',
                                      'reason': 'Contains non-carbon atoms'},
                                  {   'smiles': 'Oc1c(O)c2c3ccccc3cc3ccc4cccc1c4c23',
                                      'name': '11,12-Dihydroxybenzo[a]pyrene',
                                      'reason': 'Contains non-carbon atoms'},
                                  {   'smiles': 'Oc1ccc2ccc3c4ccccc4cc4ccc1c2c34',
                                      'name': '3-Hydroxybenzo[a]pyrene',
                                      'reason': 'Contains non-carbon atoms'},
                                  {   'smiles': 'ClC1=C(O)C2=C(O)C3=C4C(=C5C(=CC(C6=C5C=7C4=C2C(=C1O)C8=C(O)C(Cl)=C(O)C(C87)=C6O)=O)C)C(C)=CC3=O',
                                      'name': "7,7'-Dichlorohypericin",
                                      'reason': 'Contains non-carbon atoms'},
                                  {   'smiles': 'Clc1cc(Cl)c2[C@@H](c3ccccc3CCc2c1)n1ccnc1',
                                      'name': '(R)-eberconazole',
                                      'reason': 'Contains non-carbon atoms'},
                                  {   'smiles': 'Oc1cc2cccc3ccc4cc5ccccc5c1c4c23',
                                      'name': '11-Hydroxybenzo[a]pyrene',
                                      'reason': 'Contains non-carbon atoms'},
                                  {   'smiles': 'Oc1cc2cc3ccccc3c3ccc4cccc1c4c23',
                                      'name': '4-Hydroxybenzo[a]pyrene',
                                      'reason': 'Contains non-carbon atoms'}],
    'attempt': 2,
    'reasoning': 'The previous code had a few issues that led to false '
                 "positives and a low F1 score. Let's analyze the problems and "
                 'improve the code:\n'
                 '\n'
                 '1. **Problem with aromaticity detection**: The code only '
                 'checked the largest aromatic system, but some molecules have '
                 'multiple aromatic systems. We should check all aromatic '
                 'systems.\n'
                 '\n'
                 '2. **Overly permissive carbon check**: The code only checked '
                 'carbon atoms in the largest aromatic system, but non-carbon '
                 'atoms in other parts of the molecule could still make it not '
                 'a pure polycyclic arene.\n'
                 '\n'
                 '3. **Incorrect ring counting**: The code used '
                 '`CalcNumAromaticRings` which might count some non-fused '
                 'rings. We should specifically check for fused aromatic '
                 'systems.\n'
                 '\n'
                 "4. **No check for heteroatoms**: The code didn't properly "
                 'check for heteroatoms in the entire molecule.\n'
                 '\n'
                 "Here's the improved version:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 62,
    'num_false_positives': 56,
    'num_true_negatives': 142167,
    'num_false_negatives': 15,
    'num_negatives': None,
    'precision': 0.5254237288135594,
    'recall': 0.8051948051948052,
    'f1': 0.635897435897436,
    'accuracy': 0.9995010541110331,
    'negative_predictive_value': 0.9998945014136811}