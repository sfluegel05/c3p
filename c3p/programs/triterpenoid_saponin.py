"""
Classifies: CHEBI:61778 triterpenoid saponin
"""
"""
Classifies: CHEBI:75214 triterpenoid saponin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_triterpenoid_saponin(smiles: str):
    """
    Determines if a molecule is a triterpenoid saponin based on its SMILES string.
    A triterpenoid saponin is a terpene glycoside with a triterpenoid backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triterpenoid saponin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for triterpenoid backbone (typically 27-30 carbons, with multiple rings)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 27:
        return False, "Not enough carbons for a triterpenoid backbone"

    # Check for glycosidic linkages (presence of sugar moieties)
    glycosidic_pattern = Chem.MolFromSmarts("[C;!$(C=O)]-[O;!$(O=*)]-[C;!$(C=O)]")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
    if len(glycosidic_matches) == 0:
        return False, "No glycosidic linkages found"

    # Check for multiple rings (typical of triterpenoids)
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    if n_rings < 3:
        return False, "Not enough rings for a triterpenoid backbone"

    # Check for hydroxyl groups (common in triterpenoids and saponins)
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:
        return False, "Not enough hydroxyl groups for a saponin"

    # Check molecular weight (triterpenoid saponins are typically large)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for a triterpenoid saponin"

    return True, "Contains a triterpenoid backbone with glycosidic linkages"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:61778',
                          'name': 'triterpenoid saponin',
                          'definition': 'A terpene glycoside in which the '
                                        'terpene moiety is a triterpenoid.',
                          'parents': ['CHEBI:26605', 'CHEBI:61777'],
                          'xrefs': ['Wikipedia:Triterpenoid_saponin'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.5192307692307693 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'O([C@@H]1[C@@]([C@]2([C@@]([C@@]3([C@]([C@]4(C([C@]5([C@@](CC4)(CCC(C5)(C)C)C(O[C@@H]6O[C@@H]([C@H](O)[C@H](O[C@@H]7O[C@H]([C@H](O)[C@@H](O)[C@H]7O)C)[C@H]6O[C@@H]8O[C@H](C(OC9OC[C@@H](O)[C@H](O)[C@H]9O)[C@@H](O)[C@H]8O)C)C)=O)[H])=CC3)C)(CC2)C)[H])(CC1)C)[H])(C)C=O)[C@@H]%10O[C@@H]([C@@H](O)[C@H](O)[C@H]%10O[C@@H]%11O[C@@H]([C@H](O)[C@H](O)[C@H]%11O)CO)C(O)=O '
               'NAME: Goyasaponin I REASON: CORRECT Contains a triterpenoid '
               'backbone with glycosidic linkages\n'
               ' * SMILES: '
               'O([C@@H]1[C@](C2[C@@]([C@@]3([C@]([C@]4(C([C@@]5([C@@](CC4)(C[C@H](O)C(C5)(C)C)C(O)=O)[H])=CC3)C)(CC2)C)[H])(C[C@H]1O)C)(CO)C)[C@@H]6O[C@@H]([C@@H](O)[C@H](O)[C@H]6O)CO '
               'NAME: Lucyoside R REASON: CORRECT Contains a triterpenoid '
               'backbone with glycosidic linkages\n'
               ' * SMILES: '
               'C[C@@H]1O[C@@H](O[C@H]2[C@H](O)[C@@H](CO)O[C@@H](O[C@H]3[C@H](O[C@H]4CC[C@@]5(C)[C@@H](CC[C@]6(C)[C@@H]5CC=C5[C@@H]7CC(C)(C)CC[C@@]7([C@H](O)C[C@@]65C)C(=O)O[C@@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@H]5O)C4(C)C)O[C@@H]([C@@H](O)[C@@H]3O[C@@H]3OC[C@@H](O)[C@H](O)[C@H]3O)C(O)=O)[C@@H]2O)[C@H](O)[C@H](O)[C@H]1O '
               'NAME: gordonoside P REASON: CORRECT Contains a triterpenoid '
               'backbone with glycosidic linkages\n'
               ' * SMILES: '
               '[H][C@]1(O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)O[C@H]1[C@H](O)[C@@H](CO[C@@H]2O[C@H](COC(=O)C[C@@](C)(O)CC(O)=O)[C@@H](O)[C@H](O)[C@H]2O[C@]2([H])O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)O[C@@H](OC(=O)[C@]23CCC(C)(C)C[C@@]2([H])C2=CC[C@]4([H])[C@@]5(C)CC[C@H](O)[C@@](C)(C(=O)O[C@@H]6O[C@H](CO)[C@@H](O)[C@H](O)[C@H]6O)[C@]5([H])CC[C@@]4(C)[C@]2(C)CC3)[C@@H]1O '
               'NAME: dianversicoside C REASON: CORRECT Contains a '
               'triterpenoid backbone with glycosidic linkages\n'
               ' * SMILES: '
               '[H][C@@]1(CO[C@@]([H])(O[C@H]2CC[C@@]3(C)[C@@]([H])(CC[C@]4(C)[C@]3([H])CC=C3[C@]5([H])CC(C)(C)CC[C@@]5(CC[C@@]43C)C(O)=O)[C@]2(C)C=O)[C@H](O[C@]2([H])O[C@@H](C)[C@H](O)[C@@H](O)[C@H]2O)[C@H]1O)O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O '
               'NAME: salzmannianoside A REASON: CORRECT Contains a '
               'triterpenoid backbone with glycosidic linkages\n'
               ' * SMILES: '
               'O([C@@H]1[C@@]([C@]2([C@@]([C@@]3([C@]([C@]4(C([C@]5([C@@](CC4)(CC[C@](C5)(C)C(O)=O)C)[H])=CC3=O)C)(CC2)C)[H])(CC1)C)[H])(CO)C)[C@@H]6O[C@@H]([C@@H](O)[C@H](O)[C@H]6O[C@@H]7O[C@@H]([C@@H](O)[C@H](O)[C@H]7O)C(O)=O)C(O)=O '
               'NAME: Licoricesaponin g2 REASON: CORRECT Contains a '
               'triterpenoid backbone with glycosidic linkages\n'
               ' * SMILES: '
               'O(C1C(C2[C@@](C3[C@]([C@]4(C([C@@H]5[C@@](CC4)([C@@H](O[C@@H]6OC[C@H](O)[C@H](O[C@H]7O[C@@H]([C@@H](O)[C@H](O)[C@H]7O)CO)[C@H]6O)[C@@H](O)C(C5)(C)C)C)=CC3)C)(CC2)C)(CC1)C)(CO)C)[C@@H]8O[C@@H]([C@H](O)[C@H](O)[C@H]8O[C@@H]9O[C@@H]([C@H](O)[C@H](O)[C@H]9O)CO)C(O)=O '
               'NAME: Soyasaponin A2 REASON: CORRECT Contains a triterpenoid '
               'backbone with glycosidic linkages\n'
               ' * SMILES: '
               'CC1(C)CC[C@@]2(CC[C@]3(C)C(=CC[C@@H]4[C@@]5(C)CC[C@H](O[C@@H]6OC[C@H](O)[C@H](O)[C@H]6O)[C@@](C)(CO)[C@@H]5CC[C@@]34C)[C@@H]2C1)C(O)=O '
               'NAME: hederagenin 3-O-arabinoside REASON: CORRECT Contains a '
               'triterpenoid backbone with glycosidic linkages\n'
               ' * SMILES: '
               '[H][C@]1(O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)O[C@H]1[C@H](O)[C@@H](CO[C@@H]2O[C@H](COC(=O)C[C@@](C)(O)CC(O)=O)[C@@H](O)[C@H](O)[C@H]2O[C@]2([H])O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)O[C@@H](OC(=O)[C@]23CCC(C)(C)C[C@@]2([H])C2=CC[C@]4([H])[C@@]5(C)CC[C@H](O)[C@@](C)(C(O)=O)[C@]5([H])CC[C@@]4(C)[C@]2(C)CC3)[C@@H]1O '
               'NAME: dianversicoside D REASON: CORRECT Contains a '
               'triterpenoid backbone with glycosidic linkages\n'
               ' * SMILES: '
               '[H][C@]1(OC[C@@H](O[C@]2([H])O[C@H](CO)[C@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@H]1O)O[C@H]1[C@H](C)O[C@@]([H])(O[C@@H]2[C@@H](O)[C@@H](OC(=O)\\C=C/c3ccc(OC)cc3)[C@@H](C)O[C@H]2OC(=O)[C@]23CCC(C)(C)C[C@@]2([H])C2=CC[C@]4([H])[C@@]5(C)C[C@H](O)[C@H](O[C@]6([H])O[C@H](CO)[C@@H](O)[C@H](O)[C@H]6O)[C@@](C)(C(O)=O)[C@]5([H])CC[C@@]4(C)[C@]2(CO)CC3)[C@H](O)[C@@H]1O[C@]1([H])OC[C@](O)(CO)[C@H]1O '
               'NAME: Z-senegasaponin a REASON: CORRECT Contains a '
               'triterpenoid backbone with glycosidic linkages\n'
               ' * SMILES: '
               'O(C1C(C2C(C3C(C4(C(C5C(CC4)(CCC(C5)(C)C)C(O)=O)=CC3)C)(CC2)C)(CC1)C)(C)C)C6OC(C(O)C(OC7OC(C(O)C(O)C7O)CO)C6O)C(OC)=O '
               'NAME: Calenduloside G methyl ester REASON: CORRECT Contains a '
               'triterpenoid backbone with glycosidic linkages\n'
               ' * SMILES: '
               'O(C1C(C2C(C3C(C4(C(C5C(CC4)(CCC(C5)(C)C(O)=O)C)=CC3=O)C)(CC2)C)(CC1)C)(C)C)C6OC(C(O)C(O)C6O)C(O)=O '
               'NAME: Monoglucuronylglycyrrhetinic acid REASON: CORRECT '
               'Contains a triterpenoid backbone with glycosidic linkages\n'
               ' * SMILES: '
               'O[C@]1([C@]2(C=3[C@]([C@]4(C([C@@]5(C(CC4)[C@@]([C@@H](O)[C@H](O)C5)(C)C(O)=O)C)CC3)C)(CC[C@]2(CC[C@H]1C)C(O[C@@H]6O[C@@H]([C@@H](O)[C@H](O)[C@H]6O)CO)=O)C)[H])C '
               'NAME: '
               '(2r,3r,4r,6ar,6bs,8as,11r,12r,12as,14br)-2,3,12-trihydroxy-4,6a,6b,11,12,14b-hexamethyl-8a-[(2s,3r,4s,5s,6r)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxycarbonyl-1,2,3,4a,5,6,7,8,9,10,11,12a,14,14a-tetradecahydropicene-4-carboxylic '
               'acid REASON: CORRECT Contains a triterpenoid backbone with '
               'glycosidic linkages\n'
               ' * SMILES: '
               'O[C@@H]1[C@@]2(C3[C@]([C@]4(C([C@]5([C@@](CC4)(CC[C@@](C5)(C)C(O)=O)C)[H])=CC3)C)(CCC2[C@@]([C@@H](O)C1)(CO[C@@H]6O[C@H]([C@H](OC(=O)C)[C@@H](O)[C@H]6O)C)C)C)C '
               'NAME: 1alpha,3beta-Hydroxyimberbic '
               'acid-23-O-alpha-L-4-acetylrhamnopyranoside REASON: CORRECT '
               'Contains a triterpenoid backbone with glycosidic linkages\n'
               ' * SMILES: '
               '[H][C@@]1(CC[C@]2(C)[C@@H](CC[C@H]3C(C)=CC[C@@H](C(C)=C)[C@]3(C)CCC(O)=O)C(=C)CC[C@@]2([H])C1(C)C)O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1NC(C)=O '
               'NAME: lansioside A REASON: CORRECT Contains a triterpenoid '
               'backbone with glycosidic linkages\n'
               ' * SMILES: '
               'O(C1C(C2C(C3C(C4(C(C5C(C(O)C4)(CCC(C5)(C)C)C(O)=O)=CC3)C)(CC2)C)(CC1)C)(C)C)C6OC(C(O)C(O)C6O)COC7OCC(O)C(O)C7OC8OCC(O)C(O)C8O '
               'NAME: Pitheduloside D REASON: CORRECT Contains a triterpenoid '
               'backbone with glycosidic linkages\n'
               ' * SMILES: '
               '[H][C@@]12C[C@](C)(CC[C@]1(C)CC[C@]1(C)C2=CC(=O)[C@]2([H])[C@@]3(C)CC[C@H](O[C@H]4O[C@@H]([C@@H](O)[C@H](O)[C@H]4O[C@@H]4O[C@@H]([C@@H](O)[C@H](O)[C@H]4O)C(O)=O)C(O)=O)C(C)(C)[C@]3([H])CC[C@@]12C)C(O)=O '
               'NAME: glycyrrhizinic acid REASON: CORRECT Contains a '
               'triterpenoid backbone with glycosidic linkages\n'
               ' * SMILES: '
               'CC(=O)O[C@@]12CC(C)(C)CC[C@@H]1C(=O)C[C@]1(C)C2=CC[C@@H]2[C@@]3(C)CC[C@H](O[C@@H]4O[C@@H]([C@@H](O[C@@H]5O[C@H](CO)[C@H](O)[C@H](O)[C@H]5O[C@@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@H]5O)[C@H](O)[C@H]4O[C@@H]4O[C@H](CO)[C@H](O)[C@H](O)[C@H]4O)C(O)=O)C(C)(C)[C@@H]3CC[C@@]12C '
               'NAME: Camellidin I REASON: CORRECT Contains a triterpenoid '
               'backbone with glycosidic linkages\n'
               ' * SMILES: '
               'O=C(O)[C@]12[C@H](C3=CC[C@@H]4[C@@]5([C@H]([C@](C(=O)O)([C@@H](O[C@@H]6O[C@@H]([C@@H](O)[C@@H]([C@H]6O)O)CO)CC5)C)CC[C@]4([C@]3(C)CC1)C)C)C[C@](C(=O)O)(C)CC2 '
               'NAME: '
               '3beta-(beta-D-glucopyranosyloxy)olean-12-ene-23,28,30-trioic '
               'acid REASON: CORRECT Contains a triterpenoid backbone with '
               'glycosidic linkages\n'
               ' * SMILES: '
               'O(C1C(C2C(C3C(C4(C(C5C(CC4)(CCC(C5)(C)C)C(OC6OC(C(O)C(O)C6O)CO)=O)=CC3)C)(CC2)C)(CC1)C)(C)C)C7OC(C(OC8OC(C(O)C(O)C8O)CO)C(OC9OC(C(O)C(O)C9O)CO)C7O)C(O)=O '
               'NAME: '
               '6-[[4,4,6a,6b,11,11,14b-heptamethyl-8a-[3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxycarbonyl-1,2,3,4a,5,6,7,8,9,10,12,12a,14,14a-tetradecahydropicen-3-yl]oxy]-5-hydroxy-3,4-bis[[3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy]oxane-2-carboxylic '
               'acid REASON: CORRECT Contains a triterpenoid backbone with '
               'glycosidic linkages\n'
               ' * SMILES: '
               'O(C1C(C2C(C3C(C4(C(C5C(CC4)(CCC(C5)(C)C)C(O)=O)=CC3)C)(CC2)C)(CC1)C)(C)C)C6OC(C(OC7OC(C(O)C(O)C7O)CO)C(O)C6O)C(O)=O '
               'NAME: Calendulaglycoside E REASON: CORRECT Contains a '
               'triterpenoid backbone with glycosidic linkages\n'
               ' * SMILES: '
               'O([C@@H]1[C@@]([C@]2([C@@]([C@@]3([C@]([C@]4(C([C@]5([C@@](CC4)(CCC(C5)(C)C)C(O)=O)[H])=CC3)CO)(CC2)C)[H])(C[C@@H]1O)C)[H])(C)C(O)=O)[C@@H]6O[C@@H]([C@@H](O)[C@H](O)[C@H]6O)CO '
               'NAME: Tenuifolin REASON: CORRECT Contains a triterpenoid '
               'backbone with glycosidic linkages\n'
               ' * SMILES: '
               'O([C@@H]1C(C2[C@@](C3[C@]([C@]4(C([C@]5([C@@](CC4)(CCC(C5)(C)C)C(O)=O)[H])=CC3)C)(CC2)C)(CC1)C)(C)C)[C@H]6O[C@@H]([C@@H](O)[C@H](O[C@@H]7O[C@@H]([C@H](O)[C@H](O)[C@H]7O)CO)[C@H]6O)C(O)=O '
               'NAME: '
               '(2S,3S,4S,5R,6S)-6-[[(3S,6aR,6bS,8aS,12aS,14bR)-8a-carboxy-4,4,6a,6b,11,11,14b-heptamethyl-1,2,3,4a,5,6,7,8,9,10,12,12a,14,14a-tetradecahydropicen-3-yl]oxy]-3,5-dihydroxy-4-[(2S,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxane-2-carboxylic '
               'acid REASON: CORRECT Contains a triterpenoid backbone with '
               'glycosidic linkages\n'
               ' * SMILES: '
               'O([C@@H]1[C@@](C2[C@@](C3[C@]([C@]4(C([C@]5([C@@](CC4)([C@H](O)CC(C5)(C)C)C)[H])=CC3)C)(CC2)C)(CC1)C)(CO)C)[C@@H]6O[C@@H]([C@@H](O)[C@H](O)[C@H]6O[C@@H]7O[C@@H]([C@H](O)[C@H](O)[C@H]7O[C@@H]8O[C@H]([C@H](O)[C@@H](O)[C@H]8O)C)CO)C(O)=O '
               'NAME: 3-Rha(1-2)Gal(1-2)GluA-Soyasaponenol B REASON: CORRECT '
               'Contains a triterpenoid backbone with glycosidic linkages\n'
               ' * SMILES: '
               'O1[C@]2([C@]3([C@](C=4[C@]([C@]5([C@@]([C@@]6([C@@](CC5)(C([C@@H](O[C@@H]7O[C@@H]([C@@H](O)[C@H](O)[C@H]7O[C@@H]8O[C@@H]([C@@H](O)[C@H](O)[C@H]8O)C(O)=O)C(O)=O)CC6)(C)C)[H])C)(C(=O)C4)[H])C)(CC3)C)(C[C@](C2)(C1=O)C)[H])C)[H] '
               'NAME: Licoricesaponin E2 REASON: CORRECT Contains a '
               'triterpenoid backbone with glycosidic linkages\n'
               'False positives: NONE\n'
               'False negatives: SMILES: '
               'C[C@H]1[C@H]2[C@@H](C[C@@]3(C)[C@@H]4CC=C5[C@@H](C=C(O[C@@H]6O[C@H](CO)[C@@H](O)[C@H](O)[C@H]6O)C(=O)C5(C)C)[C@]4(C)C(=O)C[C@]23C)O[C@@]2(O)CC(C)(C)O[C@@]12O '
               'NAME: '
               '(16alpha,20S,22S,23S)-16,23:22,25-diepoxy-2,22,23-trihydroxycucurbita-1,5-diene-3,11-dione '
               '2-O-beta-D-glucopyranoside REASON: MISSED No carboxylic acid '
               'groups found\n'
               ' * SMILES: '
               '[H][C@]1(O[C@@H](C)[C@]([H])(O[C@@H]2OC[C@@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@H]1O)O[C@@H]1[C@@H](O)[C@@H](O)CO[C@H]1OC(=O)[C@]12CCC(C)(C)C[C@@]1([H])C1=CC[C@]3([H])[C@@]4(C)C[C@H](O)[C@H](O[C@]5([H])O[C@H](CO)[C@@H](O)[C@H](O[C@]6([H])O[C@H](CO)[C@@H](O)[C@H](O)[C@H]6O)[C@H]5O)C(CO)(CO)[C@]4([H])CC[C@@]3(C)[C@]1(C)C[C@H]2O '
               'NAME: deapio-platycodin D2 REASON: MISSED No carboxylic acid '
               'groups found\n'
               ' * SMILES: '
               'C[C@@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]2O[C@H]2C[C@]3(C)[C@H](C[C@@H](O)[C@@H]4[C@H](CC[C@@]34C)[C@@](C)(O)CCC=C(C)C)[C@@]3(C)CC[C@H](O)C(C)(C)[C@H]23)[C@H](O)[C@H](O)[C@H]1O '
               'NAME: ginsenoside Rg2 REASON: MISSED No carboxylic acid groups '
               'found\n'
               ' * SMILES: '
               '[H][C@]1(OC[C@@H](O)[C@H](O)[C@H]1O)O[C@H]1CC[C@]23C[C@]22CC[C@]4(C)[C@]([H])([C@@H](O)C[C@@]4(C)[C@]2([H])C[C@H](O[C@]2([H])O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@@]3([H])C1(C)C)[C@@]1(C)CC[C@H](O1)C(C)(C)O '
               'NAME: astragaloside IV REASON: MISSED No carboxylic acid '
               'groups found\n'
               ' * SMILES: '
               'OC1C23C4(C2)C(C5(C(CC4)(C(C(O)C5)C(CCC(O)(C(C)C)COC6OC(C(O)C(O)C6O)CO)C)C)C)CCC3C(C(O)C1)(C)C(OC7OC(C(O)C(O)C7O)CO)=O '
               'NAME: Cyclotricuspidoside C REASON: MISSED No carboxylic acid '
               'groups found\n'
               ' * SMILES: '
               'O[C@H]1[C@@]2([C@]([C@]3([C@@]([C@](CC3)([C@@H](CC[C@H](C(C)(C)O)O[C@@H]4O[C@@H]([C@@H](O)[C@H](O)[C@H]4O[C@@H]5O[C@@H]([C@@H](O)[C@H](O)[C@H]5O)CO)CO)C)[H])(C1)C)C)(CC=C6[C@]2(CC[C@H](O[C@@H]7O[C@@H]([C@@H](O[C@@H]8O[C@@H]([C@@H](O)[C@H](O)[C@H]8O)CO)[C@H](O)[C@H]7O)CO)C6(C)C)[H])[H])C '
               'NAME: isomogroside IV REASON: MISSED No carboxylic acid groups '
               'found\n'
               ' * SMILES: '
               'CO[C@H]1[C@H](O)[C@@H](CO)O[C@@H](O[C@H]2[C@H](O)[C@@H](CO)O[C@@H](O[C@@H]3[C@@H](C)O[C@@H](O[C@@H]4[C@@H](O)[C@@H](CO[C@H]4O[C@H]4CC[C@]5(C)[C@@H]6CC[C@]78[C@H](C(=C)C[C@@]7(C)C6=CC[C@H]5C4(C)C)[C@@](C)(OC8=O)C(=O)CCC(C)=C)OS(O)(=O)=O)[C@H](O)[C@H]3O)[C@@H]2O)[C@@H]1O '
               'NAME: mollisoside A REASON: MISSED No carboxylic acid groups '
               'found\n'
               ' * SMILES: '
               'C[C@H]1O[C@@H](O[C@H]2CC[C@@]3(C)[C@@H](CC[C@]4(C)[C@@H]3C=C[C@]35OC[C@@]6(CCC(C)(C)C[C@@H]36)[C@@H](O)C[C@@]45C)[C@]2(C)CO)[C@H](O)[C@@H](O[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@H]1O '
               'NAME: Saikosaponin A REASON: MISSED No carboxylic acid groups '
               'found\n'
               ' * SMILES: '
               '[H][C@@]12C[C@@H](C)[C@]3([H])[C@@]4(C)CC[C@@]56C[C@@]55CC[C@H](O[C@]7([H])OC[C@H](O)[C@H](O)[C@H]7O)C(C)(C)[C@]5([H])CC[C@@]6([H])[C@]4(C)[C@@H](O)[C@]3(O[C@@H]1C(C)(C)OC)O2 '
               'NAME: 25-O-methoxycimigenol 3-O-alpha-L-arabinopyranoside '
               'REASON: MISSED No carboxylic acid groups found\n'
               ' * SMILES: '
               'O([C@@H]1[C@@]2([C@@]([C@@]3([C@]([C@]4([C@@]([C@](CC4)([C@@](O[C@@H]5O[C@@H]([C@@H](O)[C@H](O)[C@H]5O)CO)(C/C=C\\C(O)(C)C)C)[H])([C@H](O)C3)[H])C)(C1)C)[H])(CC[C@H](O)C2(C)C)C)[H])[C@@H]6O[C@@H]([C@@H](O)[C@H](O)[C@H]6O[C@@H]7OC[C@@H](O)[C@H](O)[C@H]7O)CO '
               'NAME: Notoginsenoside H REASON: MISSED No carboxylic acid '
               'groups found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'Cl.C=1C=C(OC)C(OC)=C(C1CN2CCNCC2)OC.Cl',
                                     'name': 'Trimetazidine hydrochloride',
                                     'reason': 'Not enough carbons for a '
                                               'triterpenoid backbone'},
                                 {   'smiles': 'C1[C@@]2([C@]3(C[C@@H]([C@]4([C@]([C@@]3([C@H](C[C@@]2(C[C@@H](C1)O)[H])O)[H])(CC[C@]4([H])[C@@H](CCC(NCCS([O-])(=O)=O)=O)C)[H])C)O)[H])C',
                                     'name': 'tauroursocholate',
                                     'reason': 'Not enough carbons for a '
                                               'triterpenoid backbone'},
                                 {   'smiles': '[H]S(=O)(=O)C[C@H](N)C(O)=O',
                                     'name': 'L-cysteine-S-dioxide',
                                     'reason': 'Not enough carbons for a '
                                               'triterpenoid backbone'},
                                 {   'smiles': 'C[C@@H]1CN(C(=O)CCCN2C=C(CO[C@@H]1CN(C)CC3=CC=CC=C3OC)N=N2)[C@H](C)CO',
                                     'name': '(8R,9S)-6-[(2R)-1-hydroxypropan-2-yl]-9-[[(2-methoxyphenyl)methyl-methylamino]methyl]-8-methyl-10-oxa-1,6,13,14-tetrazabicyclo[10.2.1]pentadeca-12(15),13-dien-5-one',
                                     'reason': 'Not enough carbons for a '
                                               'triterpenoid backbone'},
                                 {   'smiles': 'Cc1ccc(cc1)N=C=Nc1ccc(C)cc1',
                                     'name': '1,3-di(p-tolyl)carbodiimide',
                                     'reason': 'Not enough carbons for a '
                                               'triterpenoid backbone'},
                                 {   'smiles': 'O=C1[C@H]([C@H]2[C@](C3=C([C@]4([C@]([C@@H]([C@@H]([C@H](O)CC(=C)C(C)C)C)CC4)(C)CC3)C)CC2)(C)CC1)C',
                                     'name': 'Gilvsin C',
                                     'reason': 'No glycosidic linkages found'},
                                 {   'smiles': 'C[C@H]1CN([C@H](COC2=C(C=C(C=C2)NC(=O)NC3=CC=CC=C3)C(=O)N(C[C@H]1OC)C)C)C(=O)C',
                                     'name': '1-[(4S,7S,8S)-5-acetyl-8-methoxy-4,7,10-trimethyl-11-oxo-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]-3-phenylurea',
                                     'reason': 'Not enough carbons for a '
                                               'triterpenoid backbone'},
                                 {   'smiles': 'C[C@@H]1CN([C@H](COC2=C(C=C(C=C2)NC(=O)C)C(=O)N(C[C@@H]1OC)C)C)C(=O)C3=NC=CN=C3',
                                     'name': 'N-[(4S,7R,8R)-8-methoxy-4,7,10-trimethyl-11-oxo-5-[oxo(2-pyrazinyl)methyl]-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]acetamide',
                                     'reason': 'Not enough carbons for a '
                                               'triterpenoid backbone'},
                                 {   'smiles': 'COc1ccc(cc1)\\C=C(NC=O)\\C(NC=O)=C\\c1ccc(OS(O)(=O)=O)cc1',
                                     'name': 'Fumiformamide',
                                     'reason': 'Not enough carbons for a '
                                               'triterpenoid backbone'},
                                 {   'smiles': 'CCCCC(CC)COCCCN',
                                     'name': '3-(2-ethylhexoxy)propan-1-amine',
                                     'reason': 'Not enough carbons for a '
                                               'triterpenoid backbone'}],
    'sample_false_negatives': [   {   'smiles': 'O=C(OC([C@@H]1O[C@](C2C3[C@@](C4(C(C5(C(C(C(OC(=O)C)CC5)(C)C)CC4)C)CC3)C)(C)CC2)(C)CC1)(C)C)C',
                                      'name': '12-Desoxydiacetylpyxinol',
                                      'reason': 'Not enough hydroxyl groups '
                                                'for a saponin'},
                                  {   'smiles': '[H][C@@]1(CC[C@]2(C)[C@]1([H])CC[C@]1([H])[C@@]3(C)CCC(=O)OC(C)(C)[C@]3([H])CC[C@@]21C)[C@]1(C)CC[C@H](O1)C(C)(C)O',
                                      'name': '20S,24S-epoxy-25-hydroxy-A-homo-4-oxadammaran-3-one',
                                      'reason': 'Not enough hydroxyl groups '
                                                'for a saponin'},
                                  {   'smiles': 'O=C(OC([C@@H]1O[C@](C2C3[C@@](C4(C(C5(C(C(C(OC(=O)C)CC5)(C)C)CC4)C)CC3O)C)(C)CC2)(C)CC1)(C)C)C',
                                      'name': 'Diacetylpyxinol',
                                      'reason': 'Not enough hydroxyl groups '
                                                'for a saponin'},
                                  {   'smiles': '[H][C@@]1(CC[C@]2(C)[C@]1([H])CC[C@]1([H])[C@@]3(C)CCC(=O)OC(C)(C)[C@]3([H])CC[C@@]21C)[C@@]1(C)CC[C@H](O1)C(C)(C)O',
                                      'name': 'rel-20R,24S-epoxy-25-hydroxy-A-homo-4-oxadammaran-3-one',
                                      'reason': 'Not enough hydroxyl groups '
                                                'for a saponin'},
                                  {   'smiles': 'O=C(OC([C@@H]1O[C@@]([C@@H]2[C@@H]3[C@@]([C@]4([C@@H]([C@@]5(C(C([C@@H](O)CC5)(C)C)CC4)C)CC3)C)(C)CC2)(C)CC1)(C)C)C',
                                      'name': '25-Acetoxy-3b-hydroxy-20(S),24(R)-epoxydammarane',
                                      'reason': 'Not enough hydroxyl groups '
                                                'for a saponin'}],
    'attempt': 1,
    'reasoning': 'The previous code had a few issues that led to false '
                 'negatives, particularly around the requirement for '
                 'carboxylic acid groups. From the error analysis, we can see '
                 'that many valid triterpenoid saponins were missed because '
                 "they didn't contain carboxylic acid groups. This suggests "
                 'that the presence of carboxylic acid groups should not be a '
                 'strict requirement for classification as a triterpenoid '
                 'saponin.\n'
                 '\n'
                 '### Key Improvements:\n'
                 '1. **Remove the Carboxylic Acid Requirement**: Many '
                 'triterpenoid saponins do not have carboxylic acid groups, so '
                 'this check should be removed.\n'
                 '2. **Refine the Glycosidic Linkage Check**: The current '
                 'pattern for glycosidic linkages might be too restrictive. We '
                 'should ensure it captures all types of glycosidic bonds.\n'
                 '3. **Adjust the Carbon Count**: While 30 carbons are typical '
                 'for triterpenoids, some structures might have slightly fewer '
                 'due to modifications.\n'
                 '4. **Improve Ring Count Check**: The ring count should be '
                 'adjusted to be more inclusive while still capturing the '
                 'triterpenoid backbone.\n'
                 '\n'
                 "Here's the improved code:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 392,
    'num_false_positives': 100,
    'num_true_negatives': 916,
    'num_false_negatives': 5,
    'num_negatives': None,
    'precision': 0.7967479674796748,
    'recall': 0.9874055415617129,
    'f1': 0.8818897637795277,
    'accuracy': 0.9256900212314225,
    'negative_predictive_value': 0.99457111834962}