"""
Classifies: CHEBI:36233 disaccharide
"""
"""
Classifies: CHEBI:36218 disaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide based on its SMILES string.
    A disaccharide is a compound in which two monosaccharides are joined by a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a disaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of at least two sugar-like rings
    # A sugar ring typically has 5 or 6 atoms, with multiple hydroxyl groups
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 1:
        return False, "No rings found, need at least one for a disaccharide"

    # Check for glycosidic bond (oxygen linking two sugar-like rings)
    glycosidic_bond_pattern = Chem.MolFromSmarts("[O;X2][C;X4][C;X4][O;X2]")
    if not mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "No glycosidic bond found between rings"

    # Count hydroxyl groups (OH) to ensure it's a sugar
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1)
    if hydroxyl_count < 4:
        return False, "Too few hydroxyl groups for a disaccharide"

    # Check molecular weight - disaccharides typically have a molecular weight between 200 and 600 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 600:
        return False, "Molecular weight outside typical range for disaccharides"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 8 or c_count > 24:
        return False, "Carbon count outside typical range for disaccharides"
    if o_count < 6 or o_count > 14:
        return False, "Oxygen count outside typical range for disaccharides"

    return True, "Contains two monosaccharide units joined by a glycosidic bond"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36233',
                          'name': 'disaccharide',
                          'definition': 'A compound in which two '
                                        'monosaccharides are joined by a '
                                        'glycosidic bond.',
                          'parents': ['CHEBI:50699'],
                          'xrefs': ['KEGG:C01911'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.7934131736526946 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'O([C@@H]1[C@H](O)[C@@H](O)C(O[C@@H]1CO)O)[C@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO '
               'NAME: alpha-D-Glcp-(1->4)-D-Galp REASON: CORRECT Contains two '
               'monosaccharide units joined by a glycosidic bond\n'
               ' * SMILES: '
               'O([C@@H]1[C@H](O)[C@@H](O)[C@@H](O[C@@H]1CO)O)[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO '
               'NAME: beta-D-Glcp-(1->4)-beta-D-Galp REASON: CORRECT Contains '
               'two monosaccharide units joined by a glycosidic bond\n'
               ' * SMILES: '
               'O([C@H]1[C@@H](O)[C@H](O)[C@H](OC1O)CO)C2O[C@@H]([C@@H](O)[C@H](O)[C@@H]2O)CO '
               'NAME: D-Manp-(1->2)-D-Manp REASON: CORRECT Contains two '
               'monosaccharide units joined by a glycosidic bond\n'
               ' * SMILES: '
               'O([C@@H]1[C@H](O)[C@@H](O)C(O[C@@H]1CO)O)[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@@H]2O)CO '
               'NAME: beta-D-Manp-(1->4)-D-Galp REASON: CORRECT Contains two '
               'monosaccharide units joined by a glycosidic bond\n'
               ' * SMILES: '
               'O([C@H]1[C@H](O)[C@@H](O)[C@H](O[C@@H]1CO)O)[C@H]2O[C@H]([C@@H](O)[C@@H](O)[C@@H]2O)C '
               'NAME: beta-L-Fucp-(1->4)-alpha-D-Glcp REASON: CORRECT Contains '
               'two monosaccharide units joined by a glycosidic bond\n'
               ' * SMILES: '
               'OC[C@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@H](O)[C@@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O '
               'NAME: alpha-D-glucosyl-(1->4)-alpha-D-mannose REASON: CORRECT '
               'Contains two monosaccharide units joined by a glycosidic bond\n'
               ' * SMILES: '
               'O([C@H]1[C@H](O)[C@H](OC(O)[C@@H]1O)CO)[C@@H]2O[C@@H]([C@H](O)[C@H](O)[C@H]2O)CO '
               'NAME: beta-D-Galp-(1->3)-D-Glcp REASON: CORRECT Contains two '
               'monosaccharide units joined by a glycosidic bond\n'
               ' * SMILES: '
               'OC[C@H]1O[C@@H](O[C@@H]2[C@@H](CO)O[C@@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O '
               'NAME: beta-cellobiose REASON: CORRECT Contains two '
               'monosaccharide units joined by a glycosidic bond\n'
               ' * SMILES: '
               'OC[C@H]1O[C@@H](O[C@@H]2[C@@H](CO)OC(O)[C@@H](O)[C@H]2O)[C@@H](O)[C@@H](O)[C@@H]1O '
               'NAME: mannobiose REASON: CORRECT Contains two monosaccharide '
               'units joined by a glycosidic bond\n'
               ' * SMILES: '
               'OC[C@H]1O[C@@H](O[C@@H]2[C@@H](CO)O[C@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O '
               'NAME: alpha-cellobiose REASON: CORRECT Contains two '
               'monosaccharide units joined by a glycosidic bond\n'
               ' * SMILES: '
               'OC[C@H]1O[C@@H](O[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O '
               'NAME: beta,beta-trehalose REASON: CORRECT Contains two '
               'monosaccharide units joined by a glycosidic bond\n'
               ' * SMILES: '
               'O([C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O)CO)[C@@H]2[C@@H](O)[C@H](O)CO[C@H]2O '
               'NAME: beta-D-Glcp-(1->2)-beta-D-Xylp REASON: CORRECT Contains '
               'two monosaccharide units joined by a glycosidic bond\n'
               ' * SMILES: '
               '[C@H]1([C@@H]([C@H]([C@@H]([C@H](O1)CO)O)O)O)O[C@H]2[C@@H]([C@@H]([C@H](O[C@@]2([C@H](CO)O)[H])O)O)O '
               'NAME: beta-D-Glcp-(1->4)-L-alpha-D-Hepp REASON: CORRECT '
               'Contains two monosaccharide units joined by a glycosidic bond\n'
               ' * SMILES: '
               'O([C@H]1[C@H](O)[C@@H](O)[C@@H](OC1)O)[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)C(O)=O '
               'NAME: beta-D-GlcpA-(1->4)-beta-D-Xylp REASON: CORRECT Contains '
               'two monosaccharide units joined by a glycosidic bond\n'
               ' * SMILES: '
               'OC[C@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@H](O)O[C@H](CO)[C@@H]2O)[C@H](O)[C@@H](O)[C@H]1O '
               'NAME: beta-(1->3)-galactobiose REASON: CORRECT Contains two '
               'monosaccharide units joined by a glycosidic bond\n'
               ' * SMILES: '
               'OC[C@H]1O[C@H](OC[C@H]2OC(O)[C@@H](O)[C@@H](O)[C@@H]2O)[C@@H](O)[C@@H](O)[C@@H]1O '
               'NAME: alpha-D-Manp-(1->6)-D-Manp REASON: CORRECT Contains two '
               'monosaccharide units joined by a glycosidic bond\n'
               ' * SMILES: '
               'O([C@H]1[C@H](O)[C@H](OC(O)[C@H]1O)CO)C2O[C@@H]([C@@H](O)[C@H](O)[C@@H]2O)CO '
               'NAME: D-Manp-(1->3)-D-Manp REASON: CORRECT Contains two '
               'monosaccharide units joined by a glycosidic bond\n'
               ' * SMILES: '
               'O1[C@@H]([C@H](O)[C@H](O)[C@@H](O)C1O)CO[C@@H]2O[C@H]([C@@H](O)[C@@H](O)[C@@H]2O)C '
               'NAME: alpha-L-Fucp-(1->6)-D-Galp REASON: CORRECT Contains two '
               'monosaccharide units joined by a glycosidic bond\n'
               ' * SMILES: '
               'O([C@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O)CO)[C@H]2[C@H](O)C(OC[C@H]2O)O '
               'NAME: alpha-D-Glcp-(1->3)-D-Arap REASON: CORRECT Contains two '
               'monosaccharide units joined by a glycosidic bond\n'
               ' * SMILES: '
               'O[C@@H]1[C@@H](O)[C@@H](O)CO[C@H]1OC[C@H]2O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]2O '
               'NAME: alpha-L-Arap-(1->6)-beta-D-Glcp REASON: CORRECT Contains '
               'two monosaccharide units joined by a glycosidic bond\n'
               ' * SMILES: '
               'O1[C@@H](OC[C@H]2OC(O)[C@H](O)[C@@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1CO '
               'NAME: beta-D-Glcp-(1->6)-D-Galp REASON: CORRECT Contains two '
               'monosaccharide units joined by a glycosidic bond\n'
               ' * SMILES: '
               'O1[C@@H](OC[C@H]2O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]2O)[C@H](O)[C@@H](O)[C@H]1CO '
               'NAME: beta-D-Xylf-(1->6)-beta-D-Glcp REASON: CORRECT Contains '
               'two monosaccharide units joined by a glycosidic bond\n'
               ' * SMILES: '
               'O([C@H]1[C@@H](O)[C@H](OC(O)[C@@H]1O)CO)[C@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO '
               'NAME: alpha-D-Glcp-(1->3)-D-Galp REASON: CORRECT Contains two '
               'monosaccharide units joined by a glycosidic bond\n'
               ' * SMILES: '
               'OC[C@H]1O[C@@H](O[C@@H]2[C@@H](CO)O[C@H](O)[C@@H](O)[C@H]2O)[C@@H](O)[C@@H](O)[C@@H]1O '
               'NAME: alpha-mannobiose REASON: CORRECT Contains two '
               'monosaccharide units joined by a glycosidic bond\n'
               ' * SMILES: '
               'O([C@@H]1[C@@H](O)[C@@H](OC1O)[C@H](O)CO)[C@@H]2O[C@H]([C@H](O)[C@H]2O)[C@H](O)CO '
               'NAME: beta-D-Galf-(1->2)-D-Galf REASON: CORRECT Contains two '
               'monosaccharide units joined by a glycosidic bond\n'
               'False positives: SMILES: '
               'CNc1ccn([C@@H]2O[C@H](COP(O)(O)=O)[C@@H](O)[C@H]2O)c(=O)n1 '
               "NAME: N(4)-methylcytidine 5'-monophosphate REASON: WRONGLY "
               'CLASSIFIED Contains two monosaccharide units joined by a '
               'glycosidic bond\n'
               ' * SMILES: O(C1CC(CCC1C)C(O)(CO)C)C2OC(C(O)C(O)C2O)CO NAME: '
               '(1S,2R,4R,8S)-p-Menthane-2,8,9-triol 2-glucoside REASON: '
               'WRONGLY CLASSIFIED Contains two monosaccharide units joined by '
               'a glycosidic bond\n'
               ' * SMILES: '
               'C=1(C=C(C=C(C1)O)CC2OC(CC2)=O)OC3[C@@H]([C@H]([C@@H]([C@H](O3)C(=O)OC)O)O)O '
               'NAME: '
               "5-(3',5'-Dihydroxyphenyl)-gamma-valerolactone-O-glucuronide-O-methyl "
               'REASON: WRONGLY CLASSIFIED Contains two monosaccharide units '
               'joined by a glycosidic bond\n'
               ' * SMILES: '
               'C12N(C(NC(C1(C)C3(C2N(C(NC3=O)=O)[C@@H]4O[C@@H]([C@H](C4)O)CO)C)=O)=O)[C@@H]5O[C@@H]([C@H](C5)O)CO '
               'NAME: thymidine dimer REASON: WRONGLY CLASSIFIED Contains two '
               'monosaccharide units joined by a glycosidic bond\n'
               ' * SMILES: '
               'COC(C)C(=O)N[C@H]1[C@H](C)OC(O[C@@H]2[C@H](O)[C@H](C)OC(O)[C@H]2O)[C@H](OC)[C@]1(C)O '
               'NAME: '
               '4,6-dideoxy-4-(2-methoxypropanamido)-3-C-Me-L-Manp2Me-(1->3)-Fucp '
               'REASON: WRONGLY CLASSIFIED Contains two monosaccharide units '
               'joined by a glycosidic bond\n'
               ' * SMILES: '
               'O1C(C(O)C2=C1C(OC3OC(C(O)C(O)C3O)C(O)=O)=C4OC(=O)C=CC4=C2)C(O)(C)C '
               'NAME: '
               '3,4,5-trihydroxy-6-{[3-hydroxy-2-(2-hydroxypropan-2-yl)-7-oxo-2H,3H,7H-furo[3,2-g]chromen-9-yl]oxy}oxane-2-carboxylic '
               'acid REASON: WRONGLY CLASSIFIED Contains two monosaccharide '
               'units joined by a glycosidic bond\n'
               ' * SMILES: O1C(C(O)C(O)C(O)C1OC(=O)C=C(C)C)COC2OCC(O)(C2O)CO '
               'NAME: 1-(3-Methyl-2-butenoyl)-6-apiosylglucose REASON: WRONGLY '
               'CLASSIFIED Contains two monosaccharide units joined by a '
               'glycosidic bond\n'
               ' * SMILES: '
               'C123C45[C@@H](OC([C@@H]4O)=O)O[C@]1([C@]6([C@@H](C(=O)O[C@]6([C@@H]2O)[H])C)O)C(=O)O[C@@]3([C@H](O)[C@]5(C(C)(C)C)[H])[H] '
               'NAME: Ginkgolide C REASON: WRONGLY CLASSIFIED Contains two '
               'monosaccharide units joined by a glycosidic bond\n'
               ' * SMILES: '
               'O1C(OC=2C=C3[O+]=C(C(O)=CC3=C(O)C2)C4=CC(O)=C(O)C=C4)C(O)C(O)C(O)C1 '
               'NAME: Cyanidin 7-arabinoside REASON: WRONGLY CLASSIFIED '
               'Contains two monosaccharide units joined by a glycosidic bond\n'
               ' * SMILES: '
               'O1[C@@H](OC2=C(O)C=C([C@@H](O)CNC)C=C2)[C@H](O)[C@@H](O)[C@H](O)[C@H]1C(O)=O '
               'NAME: Epinephrine glucuronide REASON: WRONGLY CLASSIFIED '
               'Contains two monosaccharide units joined by a glycosidic bond\n'
               ' * SMILES: '
               'O1C(C(O)C(O)C(O)C1OC2=CC=C(C=C2)C=CC3=CC(O)=CC(O)=C3)C(O)=O '
               "NAME: Resveratrol 4'-O-D-glucuronide REASON: WRONGLY "
               'CLASSIFIED Contains two monosaccharide units joined by a '
               'glycosidic bond\n'
               ' * SMILES: '
               'Nc1nc2NCC(CNc3ccc(cc3)[C@@H]3O[C@H](COP([O-])([O-])=O)[C@@H](O)[C@H]3O)=Nc2c(=O)[nH]1 '
               'NAME: '
               'N-[(7,8-dihydropterin-6-yl)methyl]-4-(beta-D-ribofuranosyl)aniline '
               "5'-phosphate(2-) REASON: WRONGLY CLASSIFIED Contains two "
               'monosaccharide units joined by a glycosidic bond\n'
               ' * SMILES: O1C(C(O)C(O)C(O)C1OC2=C(OC)C=C(C=C2)C=C)C(O)=O '
               'NAME: '
               '6-(4-ethenyl-2-methoxyphenoxy)-3,4,5-trihydroxyoxane-2-carboxylic '
               'acid REASON: WRONGLY CLASSIFIED Contains two monosaccharide '
               'units joined by a glycosidic bond\n'
               ' * SMILES: '
               'C[C@H](O)\\C=C\\[C@@]1(O)C(C)(C)CCC[C@@]1(C)O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O '
               'NAME: Rehmaionoside A REASON: WRONGLY CLASSIFIED Contains two '
               'monosaccharide units joined by a glycosidic bond\n'
               ' * SMILES: '
               '[Na+].C[C@@H]1O[C@H](O)[C@H](O[C@H]2OC(=C[C@H](O)[C@H]2O)C([O-])=O)[C@H](O)[C@H]1O '
               'NAME: Sodium '
               '2-O-L-rhamnopyranosyl-4-deoxy-alpha-L-threo-hex-4-eno-pyranosiduronate '
               'REASON: WRONGLY CLASSIFIED Contains two monosaccharide units '
               'joined by a glycosidic bond\n'
               ' * SMILES: '
               'COC1=CC=CC=C1S(=O)(=O)N2C[C@H](COC[C@H]3[C@H]2CC[C@H](O3)CC(=O)OC)O '
               'NAME: '
               '2-[(3R,6aR,8S,10aR)-3-hydroxy-1-(2-methoxyphenyl)sulfonyl-3,4,6,6a,8,9,10,10a-octahydro-2H-pyrano[2,3-c][1,5]oxazocin-8-yl]acetic '
               'acid methyl ester REASON: WRONGLY CLASSIFIED Contains two '
               'monosaccharide units joined by a glycosidic bond\n'
               ' * SMILES: '
               'CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@H](O)[C@@H]1O[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O '
               'NAME: beta-D-Galp-(1->3)-beta-D-GalpNAc REASON: WRONGLY '
               'CLASSIFIED Contains two monosaccharide units joined by a '
               'glycosidic bond\n'
               ' * SMILES: '
               'CNCc1cc(O)cc(O[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)c1O '
               'NAME: zinolol REASON: WRONGLY CLASSIFIED Contains two '
               'monosaccharide units joined by a glycosidic bond\n'
               ' * SMILES: '
               'O=C\\1OCC(/C1=C(\\O[C@@H]2O[C@H]([C@H](O)[C@H]([C@H]2O)O)C)/C=C/C)CO '
               'NAME: Rhamnosyllactone B1 REASON: WRONGLY CLASSIFIED Contains '
               'two monosaccharide units joined by a glycosidic bond\n'
               ' * SMILES: '
               'CC(=O)N[C@H]1C(O)O[C@H](CO)[C@H](O)[C@@H]1O[C@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1NC(C)=O '
               'NAME: alpha-D-GalpNAc-(1->3)-D-GalpNAc REASON: WRONGLY '
               'CLASSIFIED Contains two monosaccharide units joined by a '
               'glycosidic bond\n'
               ' * SMILES: '
               'O1[C@@H]([C@@H](O)[C@H](O)[C@@H](O)[C@H]1OC2=CC=3OC(=O)C=C(C3C=C2)C)CO '
               'NAME: 4-Methylumbelliferyl-alpha-D-glucopyranoside REASON: '
               'WRONGLY CLASSIFIED Contains two monosaccharide units joined by '
               'a glycosidic bond\n'
               ' * SMILES: '
               'C\\C=C(\\C)C(=O)O[C@@H]1C[C@]2(CO2)[C@@H]2[C@H](O)[C@H]3O[C@@]3(C)[C@]2(O)[C@H]2OC(=O)C(=C)[C@H]12 '
               'NAME: Eupatoroxin REASON: WRONGLY CLASSIFIED Contains two '
               'monosaccharide units joined by a glycosidic bond\n'
               ' * SMILES: '
               'CNCC(=O)N[C@H](CO)C(=O)N[C@H]1[C@H](O)[C@@H](O)[C@@H](O[C@@H]1C(N)=O)N1C=CC(N)=NC1=O '
               'NAME: gougerotin REASON: WRONGLY CLASSIFIED Contains two '
               'monosaccharide units joined by a glycosidic bond\n'
               ' * SMILES: '
               'O([C@@H]1[C@@H](NC(=O)C)[C@H](O[C@@H]([C@@H]1O)CO)O)C2O[C@@H]([C@@H](O)[C@H](O)[C@H]2NC(=O)C)CO '
               'NAME: '
               'N-[(3R,4R,5S,6R)-2-[(2S,3R,4R,5R,6R)-3-Acetamido-2,5-dihydroxy-6-(hydroxymethyl)oxan-4-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]acetamide '
               'REASON: WRONGLY CLASSIFIED Contains two monosaccharide units '
               'joined by a glycosidic bond\n'
               ' * SMILES: '
               'S1C(N)(C(=O)NC=2C(=O)OC=3C(OC)=C(OC)C=CC3C2)C[C@@]4([C@H]1C=C[C@@H](O)[C@@H]4O)O '
               'NAME: Aspergillazine C REASON: WRONGLY CLASSIFIED Contains two '
               'monosaccharide units joined by a glycosidic bond\n'
               'False negatives: SMILES: '
               'OC[C@@H](O)[C@@H](O[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)[C@H](O)C=O '
               'NAME: 3-O-beta-D-galactopyranosyl-D-arabinose REASON: MISSED '
               'Less than 2 rings found, need exactly 2 for a disaccharide\n'
               ' * SMILES: '
               '[H][C@@]1(O[C@](O)(C[C@@H](O[C@]2(O[C@]([H])([C@H](O)CO)[C@H](O)[C@H](O)[C@@H]2O)C(O)=O)[C@H]1O)C(O)=O)[C@H](O)CO '
               'NAME: alpha-D-Ko-(2->4)-alpha-D-Kdo REASON: MISSED Oxygen '
               'count outside typical range for disaccharides\n'
               ' * SMILES: '
               'O1[C@@H]([C@@H](O[C@H]2OC(=C(O)[C@H](O)[C@H]2OS(O)(=O)=O)C(O)=O)[C@H](O)[C@@H](NS(O)(=O)=O)[C@@H]1O)CO '
               'NAME: '
               '2-O-sulfo-beta-L-threo-hex-4-enopyranuronosyl-(1->4)-beta-D-GlcpNS '
               'REASON: MISSED Molecular weight outside typical range for '
               'disaccharides\n'
               ' * SMILES: '
               '[H][C@]1(CC=C[C@@](OC[C@@H](O)[C@@]2([H])O[C@](O)(C[C@@H](O)[C@H]2O)C(O)=O)(O1)C(O)=O)[C@H](O)CO '
               'NAME: alpha-3,4-dehydro-3,4,5-trideoxy-Kdo-(2->8)-alpha-Kdo '
               'REASON: MISSED Oxygen count outside typical range for '
               'disaccharides\n'
               ' * SMILES: '
               'S(O)(=O)(=O)N[C@@H]1[C@@H](O)[C@H](O[C@@H]2OC(=C(O)[C@H](O)[C@H]2O)C(O)=O)[C@H](OC1O)CO '
               'NAME: '
               '(2R,3R,4R)-2-[(2R,3S,4R,5R)-4,6-Dihydroxy-2-(hydroxymethyl)-5-(sulfoamino)oxan-3-yl]oxy-3,4,5-trihydroxy-3,4-dihydro-2H-pyran-6-carboxylic '
               'acid REASON: MISSED Oxygen count outside typical range for '
               'disaccharides\n'
               ' * SMILES: '
               'OC[C@H]1O[C@@H](O[C@H]([C@@H](O)C=O)[C@H](O)[C@H](O)C(O)=O)[C@H](O)[C@@H](O)[C@@H]1O '
               'NAME: 3-O-beta-D-glucosyl-D-glucuronic acid REASON: MISSED '
               'Less than 2 rings found, need exactly 2 for a disaccharide\n'
               ' * SMILES: '
               'OC[C@@H](O)[C@@H](O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)[C@H](O)C=O '
               'NAME: beta-D-glucosyl-(1->4)-aldehydo-D-mannose REASON: MISSED '
               'Less than 2 rings found, need exactly 2 for a disaccharide\n'
               ' * SMILES: '
               'C([C@@H]1[C@@H]([C@@H]([C@H]([C@H](O1)OC2[C@@H]([C@H](C([C@H]([C@H]2O)O)O)O)O)O)O)O)O.O.O '
               'NAME: Galactinol dihydrate REASON: MISSED Oxygen count outside '
               'typical range for disaccharides\n'
               ' * SMILES: '
               'O1[C@]2([C@@H](O)[C@H](O)[C@@H](O)[C@@]1(OCC[C@H](CCC=C(C(OC[C@]3(O)[C@@H](O)[C@@](OC3)(OC2)[H])=O)C)C)[H])[H] '
               'NAME: Urceolide REASON: MISSED Carbon count outside typical '
               'range for disaccharides\n'
               ' * SMILES: '
               'S(O)(=O)(=O)N[C@@H]1[C@@H](O)[C@H](O[C@@H]2OC(=C(O)[C@H](O)[C@H]2OS(O)(=O)=O)C(O)=O)[C@H](OC1O)CO '
               'NAME: '
               '(2R,3R,4R)-2-[(2R,3S,4R,5R)-4,6-Dihydroxy-2-(hydroxymethyl)-5-(sulfoamino)oxan-3-yl]oxy-4,5-dihydroxy-3-sulfooxy-3,4-dihydro-2H-pyran-6-carboxylic '
               'acid REASON: MISSED Molecular weight outside typical range for '
               'disaccharides\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'O(C1C(O)C(OC(OC2=C(OC=3C(C2=O)=C(O)C=C(O)C3)C4=CC(O)=C(O)C=C4)C1O)CO)C5OC(C(O)C(O)C5O)CO',
                                     'name': 'Quercetin 3-beta-laminaribioside',
                                     'reason': 'Molecular weight outside '
                                               'typical range for '
                                               'disaccharides'},
                                 {   'smiles': 'O=C(N[C@@H](CC1=CC=C(O)C=C1)C(O)=O)[C@H]2N(CCC2)C(=O)[C@@H](N)CC3=CC=CC=C3',
                                     'name': 'Phe-Pro-Tyr',
                                     'reason': 'No glycosidic bond found '
                                               'between rings'},
                                 {   'smiles': 'O[C@@H]1[C@@H](C(=C2C[C@](CO)(C)C[C@H]2[C@H](C1)C)CO)CO',
                                     'name': '4alpha,11,12,14-tetrahydroxy-1-tremulene',
                                     'reason': 'No glycosidic bond found '
                                               'between rings'},
                                 {   'smiles': 'O(CCCCCCCCCC\\C=C/C=C\\CC)C(=O)C',
                                     'name': '11Z,13Z-Hexadecadienyl acetate',
                                     'reason': 'No rings found, need at least '
                                               'one for a disaccharide'},
                                 {   'smiles': 'O=C(NC(C(=O)NC(CN(CCO)C)C)(C)C)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C1N(C(=O)C(CCCCCCCC)C)CCC1)CC(CC(O)CC(=O)CC)C)C)(C)C)C(CC)C)C(C)C)(C)C',
                                     'name': 'Roseoferin A1',
                                     'reason': 'No glycosidic bond found '
                                               'between rings'},
                                 {   'smiles': 'O=C1C(=C(CC2=CC=C(O)C=C2)C(C1(C3=CC=C(O)C=C3)C4(C5=CC=C(O)C=C5)C(=O)C(CC6=CC=C(O)C=C6)=C(C4=O)CC7=CC=C(O)C=C7)=O)CC8=CC=C(O)C=C8',
                                     'name': 'Nostotrebin 6',
                                     'reason': 'No glycosidic bond found '
                                               'between rings'},
                                 {   'smiles': 'ClC1=C(O)C=CC(=C1)NC2=NC=NC=3C2=CC(OCCCNCCCO)=C(C3)OC',
                                     'name': '2-chloro-4-[(6-{3-[(3-hydroxypropyl)amino]propoxy}-7-methoxyquinazolin-4-yl)amino]phenol',
                                     'reason': 'No glycosidic bond found '
                                               'between rings'},
                                 {   'smiles': 'C1[C@@H]2C([C@@H](N2)CN1CC3=COC=N3)C4=CC=C(C=C4)C5=CC=C(C=C5)C#N',
                                     'name': '4-[4-[(1S,5R)-3-(4-oxazolylmethyl)-3,6-diazabicyclo[3.1.1]heptan-7-yl]phenyl]benzonitrile',
                                     'reason': 'No glycosidic bond found '
                                               'between rings'},
                                 {   'smiles': 'S(CC[C@H](NC(=O)[C@@H](NC(=O)[C@@H](N)C)[C@H](CC)C)C(O)=O)C',
                                     'name': 'Ala-Ile-Met',
                                     'reason': 'No rings found, need at least '
                                               'one for a disaccharide'},
                                 {   'smiles': 'O=C([C@H]1C(=C[C@@H](O)[C@@H]2[C@@H]1CC[C@@H](C2)C)C)CCOC(=O)C',
                                     'name': 'Pallidopenilline G',
                                     'reason': 'No glycosidic bond found '
                                               'between rings'}],
    'sample_false_negatives': [   {   'smiles': '[H][C@@]1(O[C@](O)(C[C@@H](O[C@]2(O[C@]([H])([C@H](O)CO)[C@H](O)[C@H](O)[C@@H]2O)C(O)=O)[C@H]1O)C(O)=O)[C@H](O)CO',
                                      'name': 'alpha-D-Ko-(2->4)-alpha-D-Kdo',
                                      'reason': 'Oxygen count outside typical '
                                                'range for disaccharides'},
                                  {   'smiles': 'O1[C@@H]([C@@H](O[C@H]2OC(=C(O)[C@H](O)[C@H]2OS(O)(=O)=O)C(O)=O)[C@H](O)[C@@H](NS(O)(=O)=O)[C@@H]1O)CO',
                                      'name': '2-O-sulfo-beta-L-threo-hex-4-enopyranuronosyl-(1->4)-beta-D-GlcpNS',
                                      'reason': 'Oxygen count outside typical '
                                                'range for disaccharides'},
                                  {   'smiles': 'S(O)(=O)(=O)N[C@@H]1[C@@H](O)[C@H](O[C@@H]2OC(=C(O)[C@H](O)[C@H]2OS(O)(=O)=O)C(O)=O)[C@H](OC1O)CO',
                                      'name': '(2R,3R,4R)-2-[(2R,3S,4R,5R)-4,6-Dihydroxy-2-(hydroxymethyl)-5-(sulfoamino)oxan-3-yl]oxy-4,5-dihydroxy-3-sulfooxy-3,4-dihydro-2H-pyran-6-carboxylic '
                                              'acid',
                                      'reason': 'Oxygen count outside typical '
                                                'range for disaccharides'},
                                  {   'smiles': 'C[C@@H](CC(O)=O)CC(=O)OC[C@H](C)\\C=C(/C)C[C@@H]1O[C@]2(C)C[C@@H]3O[C@@H]4C[C@@H]5O[C@@H]6CC[C@@H]7O[C@@H]8C[C@@H]9O[C@H](C[C@H](O)C[C@@H]%10C[C@@H](C)[C@H](O%10)[C@@H](C)CC(O)=O)[C@@H](O)CC[C@H]9O[C@@]8(C)C[C@@]7(C)O[C@H]6C\\C=C/[C@@H](C)[C@H]5O[C@@]4(C)[C@H](O)[C@H]3O[C@H]2CC1=C',
                                      'name': 'Gambieric acid C',
                                      'reason': 'Molecular weight outside '
                                                'typical range for '
                                                'disaccharides'},
                                  {   'smiles': '[H][C@]1(OC[C@H](O)[C@H]1OC(=O)c1cc(O)cc(O)c1)[C@]1(CC[C@H](Cc2ccc(O)c(O)c2)C(=O)O1)c1ccc(O)c(O)c1',
                                      'name': 'chimarrhinin, (rel)-',
                                      'reason': 'Carbon count outside typical '
                                                'range for disaccharides'},
                                  {   'smiles': 'CC(C)(O)[C@H]1C[C@@]2(O)C=C(CC=C)C(=O)C(CC3=C4O[C@H](C[C@@]4(O)C=C(CC=C)C3=O)C(C)(C)O)=C2O1',
                                      'name': 'Illicidione B',
                                      'reason': 'Carbon count outside typical '
                                                'range for disaccharides'},
                                  {   'smiles': 'CC(C)(O)[C@H]1C[C@@]2(O)C[C@H](CC=C)C(=O)C(CC3=C4O[C@H](C[C@@]4(O)C[C@H](CC=C)C3=O)C(C)(C)O)=C2O1',
                                      'name': 'Illicidione C',
                                      'reason': 'Carbon count outside typical '
                                                'range for disaccharides'},
                                  {   'smiles': '[H][C@@]1(O[C@](O)([C@@H](O)[C@@H](O[C@@]2(C[C@@H](O)[C@@H](O)[C@]([H])(O2)[C@H](O)CO)C(O)=O)[C@H]1O)C(O)=O)[C@H](O)CO',
                                      'name': 'alpha-D-Kdo-(2->4)-alpha-D-Ko',
                                      'reason': 'Oxygen count outside typical '
                                                'range for disaccharides'},
                                  {   'smiles': 'S(O)(=O)(=O)N[C@@H]1[C@@H](O)[C@H](O[C@@H]2OC(=C(O)[C@H](O)[C@H]2OS(O)(=O)=O)C(O)=O)[C@H](O[C@@H]1O)CO',
                                      'name': '(2R,3R,4R)-2-[(2R,3S,4R,5R,6S)-4,6-Dihydroxy-2-(hydroxymethyl)-5-(sulfoamino)oxan-3-yl]oxy-4,5-dihydroxy-3-sulfooxy-3,4-dihydro-2H-pyran-6-carboxylic '
                                              'acid',
                                      'reason': 'Oxygen count outside typical '
                                                'range for disaccharides'},
                                  {   'smiles': 'C[C@@H](CC(O)=O)CC(=O)OC[C@H](C)\\C=C(/C)C[C@@H]1O[C@]2(C)C[C@@H]3O[C@@H]4C[C@@H]5O[C@@H]6CC[C@@H]7O[C@@H]8C[C@@H]9O[C@H](C[C@H](O)C[C@@H]%10C[C@@H](C)[C@H](O%10)[C@@H](C)CC(O)=O)[C@@](C)(O)CC[C@H]9O[C@@]8(C)C[C@@]7(C)O[C@H]6C\\C=C/[C@@H](C)[C@H]5O[C@@]4(C)[C@H](O)[C@H]3O[C@H]2CC1=C',
                                      'name': 'Gambieric acid D',
                                      'reason': 'Molecular weight outside '
                                                'typical range for '
                                                'disaccharides'}],
    'attempt': 1,
    'reasoning': 'The previous code had a few issues that led to false '
                 'positives and false negatives:\n'
                 '\n'
                 '1. **False Positives**: The code incorrectly classified some '
                 'molecules as disaccharides because:\n'
                 '   - The glycosidic bond pattern `[O;X2][C;X4][C;X4][O;X2]` '
                 'was too generic and matched non-glycosidic bonds.\n'
                 '   - The molecular weight and oxygen/carbon count ranges '
                 'were too broad, allowing non-disaccharides to pass.\n'
                 '\n'
                 '2. **False Negatives**: The code missed some true '
                 'disaccharides because:\n'
                 '   - The ring count check was too strict. Some disaccharides '
                 'may have additional rings or may not have exactly 2 rings.\n'
                 '   - The molecular weight and oxygen/carbon count ranges '
                 'were too restrictive, excluding valid disaccharides.\n'
                 '\n'
                 'To improve the code, we can:\n'
                 '- Use a more specific glycosidic bond pattern that matches '
                 'the typical structure of a glycosidic bond between two sugar '
                 'rings.\n'
                 '- Relax the ring count check to allow for additional rings '
                 'or slightly different ring structures.\n'
                 '- Adjust the molecular weight and oxygen/carbon count ranges '
                 'to be more inclusive of valid disaccharides.\n'
                 '- Add additional checks for the presence of sugar-like '
                 'structures (e.g., multiple hydroxyl groups, specific ring '
                 'sizes).\n'
                 '\n'
                 "Here's the improved code:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 291,
    'num_false_positives': 100,
    'num_true_negatives': 1931,
    'num_false_negatives': 12,
    'num_negatives': None,
    'precision': 0.7442455242966752,
    'recall': 0.9603960396039604,
    'f1': 0.8386167146974064,
    'accuracy': 0.9520137103684662,
    'negative_predictive_value': 0.9938239835306227}