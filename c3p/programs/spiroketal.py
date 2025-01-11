"""
Classifies: CHEBI:72600 spiroketal
"""
from rdkit import Chem

def is_spiroketal(smiles: str):
    """
    Determines if a molecule is a spiroketal based on its SMILES string.
    A spiroketal is defined as a spiro center where a carbon is the only shared atom between two rings,
    each containing at least one oxygen atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a spiroketal, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for spiro center
    # Look for a carbon spiro center that is connected to two distinct rings
    # Each ring should contain an oxygen atom
    spiroketal_smarts = "[C]([O])([*])"
    spiroketal_pattern = Chem.MolFromSmarts(spiroketal_smarts)
    if not spiroketal_pattern:
        return None, "Invalid SMARTS pattern"

    # Check for pattern in the molecule
    matches = mol.GetSubstructMatches(spiroketal_pattern)
    ring_info = mol.GetRingInfo()

    for match in matches:
        carbon_idx, o1_idx, _ = match
        
        # Ensure the carbon atom is part of exactly two rings
        rings_with_carbon = [ring for ring in ring_info.AtomRings() if carbon_idx in ring]
        if len(rings_with_carbon) != 2:
            continue
        
        # Verify that the oxygen atom is part of only one of the two rings
        o1_rings = [ring for ring in rings_with_carbon if o1_idx in ring]
        if len(o1_rings) != 1:
            continue
        
        return True, "Found a spiroketal with a carbon connecting two distinct ring oxygens"

    return False, "No spiroketal structure found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:72600',
                          'name': 'spiroketal',
                          'definition': 'A cyclic ketal in which the ketal '
                                        'carbon is the only common atom of two '
                                        'rings.',
                          'parents': ['CHEBI:37948', 'CHEBI:59779'],
                          'xrefs': [   'PMID:19262920',
                                       'PMID:20024126',
                                       'PMID:21076755',
                                       'PMID:21604735',
                                       'PMID:21860857',
                                       'PMID:22421755'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.6555183946488294 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'O=C1O[C@@H]2C[C@]3(O[C@@H](C2)CC=C(C[C@H](C=CC=C([C@]4([C@H]1C=C([C@@H](O)[C@H]4O)C)O)C)C)C)O[C@@H]([C@@H](C)CC3)C '
               'NAME: Milbemycin beta11 REASON: CORRECT Found a spiroketal '
               'with a carbon connecting two distinct ring oxygens\n'
               ' * SMILES: '
               'O=C1O[C@@H]2C[C@@]3(O[C@H](/C(=C/C)/C)[C@@H](C)C[C@H]3O)O[C@@H](C2)CC=C(C[C@H](C=CC=C([C@]4([C@H]1C=C([C@@H](OC)C4O)C)O)C)C)C '
               'NAME: VM44867 REASON: CORRECT Found a spiroketal with a carbon '
               'connecting two distinct ring oxygens\n'
               ' * SMILES: '
               'CC[C@H]([C@H]1CC[C@H](C)[C@@H](O1)[C@@H](C)[C@H](O)[C@H](C)C(=O)[C@H](CC)[C@H]1O[C@@]2(O[C@@]3(CC[C@](C)(O3)C3CC[C@](O)(CC)[C@H](C)O3)[C@H](O)C=C2)[C@H](C)C[C@@H]1C)C(O)=O '
               'NAME: Salinomycin REASON: CORRECT Found a spiroketal with a '
               'carbon connecting two distinct ring oxygens\n'
               ' * SMILES: '
               'O1C23O[C@](C[C@]1(CC2)C=O)([C@]4(O[C@@]([C@@H](O)[C@]5(OC6(O[C@@](CCC6)([C@H](C(OC7C[C@](O[C@@]7(C=CC(=C[C@@H](C[C@@]8(O[C@@]3(CC8)[H])C)C)C)[H])([C@]9(OCC[C@H]([C@H]9O)C)O)[H])=O)C)[H])CC5)[H])(CC4=O)C)[H])[H] '
               'NAME: Pectenotoxin 3 REASON: CORRECT Found a spiroketal with a '
               'carbon connecting two distinct ring oxygens\n'
               ' * SMILES: '
               'O=C1O[C@@H]2C[C@]3(O[C@@H](C2)CC=C(C[C@H](C=CC=C4[C@]5([C@H]1C=C([C@@H](O)[C@H]5OC4)C)O)C)C)O[C@@H]([C@@H](C)[C@@H]([C@H]3O)OC(=O)C(CC)C)CC '
               'NAME: Milbemycin alpha30 REASON: CORRECT Found a spiroketal '
               'with a carbon connecting two distinct ring oxygens\n'
               ' * SMILES: '
               'O1C23O[C@](C[C@]1(CC2)C(O)=O)([C@]4(O[C@@]([C@@H](O)[C@]5(OC6(O[C@@](CCC6)([C@H](C(O[C@@]7(C[C@](O[C@@]7(C=CC(=C[C@@H](C[C@@]8(O[C@@]3(CC8)[H])C)C)C)[H])([C@]9(OCC[C@H]([C@H]9O)C)O)[H])[H])=O)C)[H])CC5)[H])(CC4=O)C)[H])[H] '
               'NAME: Pectenotoxin 7 REASON: CORRECT Found a spiroketal with a '
               'carbon connecting two distinct ring oxygens\n'
               ' * SMILES: '
               '[C@@]12(C3=CC=C[C@@H](CC(=CC[C@]4(O[C@@]5(C[C@H](C4)OC([C@]1([H])C=C(C)[C@H]([C@@]2([H])OC3)O)=O)O[C@@H]([C@H](C[C@H]5O)C)/C(/C)=C/C)[H])C)C)O '
               'NAME: VM 44866 REASON: CORRECT Found a spiroketal with a '
               'carbon connecting two distinct ring oxygens\n'
               ' * SMILES: '
               'O=C1O[C@@H]2C[C@@]3(O[C@H](/C(=C/C)/C)[C@@H](C)C[C@H]3O)O[C@@H](C2)CC=C([C@H](O)[C@H](C=CC=C4[C@]5([C@H]1C=C([C@@H](OC)[C@H]5OC4)C)O)C)C '
               'NAME: UK-78618 REASON: CORRECT Found a spiroketal with a '
               'carbon connecting two distinct ring oxygens\n'
               ' * SMILES: '
               '[C@@]12(C3=CC=C[C@@H](CC(=CC[C@]4(O[C@@]5(C[C@H](C4)OC([C@]1([H])C=C(COC(/C(=C/C)/C)=O)[C@H]([C@@]2([H])OC3)OC)=O)O[C@@H]([C@H](C[C@H]5O)C)/C(/C)=C/C)[H])C)C)O '
               'NAME: VM 44865 REASON: CORRECT Found a spiroketal with a '
               'carbon connecting two distinct ring oxygens\n'
               ' * SMILES: '
               'O=C1OC2CC3(OC(C(C)CC3)C(CC)C)OC(C2)CC=C(C(OC4OC(C(OC5OC(C(O)C(C5)OC)C)C(C4)OC)C)C(C=CC=C6C7(C1C=C(C(O)C7OC6)C)O)C)C '
               'NAME: epi-Ivermectin B1a REASON: CORRECT Found a spiroketal '
               'with a carbon connecting two distinct ring oxygens\n'
               ' * SMILES: '
               'O=C1O[C@@H]2C[C@]3(O[C@@H](C2)CC=C(C[C@H](C=CC=C4[C@]5([C@H]1C=C([C@@H](O)[C@H]5OC4)COC(=O)C=6NC=CC6)O)C)C)O[C@@H]([C@@H](C)[C@@H]([C@H]3O)OC(=O)C(CC(C)C)C)C '
               'NAME: Milbemycin alpha29 REASON: CORRECT Found a spiroketal '
               'with a carbon connecting two distinct ring oxygens\n'
               ' * SMILES: '
               'O=C1O[C@@H]2C[C@@]3(O[C@H](C(=CC)C)[C@@H](C)C[C@H]3O)O[C@@H](C2)CC=C(C[C@H](C=CC=C4[C@]5([C@H]1C=C([C@@H](OC)[C@H]5OC4)C)O)C)C '
               'NAME: VM-44864 REASON: CORRECT Found a spiroketal with a '
               'carbon connecting two distinct ring oxygens\n'
               ' * SMILES: '
               'CC[C@H](C)CC\\C=C\\C=C(/C)[C@@H](O)C\\C=C\\C=C\\C(=O)O[C@@H]1[C@@H](O)[C@@]2(OCc3cc(O)cc(O)c23)O[C@H](CO)[C@H]1O[C@@H]1O[C@H](COC(=O)\\C=C\\C=C/C=C/C(O)CC)[C@H](O)[C@H](O)[C@H]1O '
               'NAME: papulacandin B REASON: CORRECT Found a spiroketal with a '
               'carbon connecting two distinct ring oxygens\n'
               ' * SMILES: '
               'O=C1O[C@@H]2C[C@]3(O[C@@H](C2)CC=C(C[C@H](C=CC=C([C@]4([C@H]1C=C([C@@H](OC)[C@H]4O)C)O)CO)C)C)O[C@@H]([C@@H](C)CC3)C '
               'NAME: Milbemycin beta9 REASON: CORRECT Found a spiroketal with '
               'a carbon connecting two distinct ring oxygens\n'
               ' * SMILES: '
               'O=C1OC2CC3(OC(C(=CC(C)C)C)C(C)C(C3)O)OC(C2)CC=C(CC(C=CC=C4C5(C1C=C(C(O)C5OC4)C)O)C)C '
               'NAME: Nemadectin REASON: CORRECT Found a spiroketal with a '
               'carbon connecting two distinct ring oxygens\n'
               ' * SMILES: '
               'O=C1O[C@@H]2C[C@@]3(O[C@H](/C(=C/C)/C)[C@@H](C)[C@@H]([C@H]3O)OC(=O)C(C)C)O[C@@H](C2)CC=C([C@H](O)[C@H](C=CC=C4[C@]5([C@H]1C=C([C@@H](OC)[C@H]5OC4)C)O)C)C '
               'NAME: UK-78614 REASON: CORRECT Found a spiroketal with a '
               'carbon connecting two distinct ring oxygens\n'
               ' * SMILES: '
               'O=C1O[C@@H]2C[C@@]3(O[C@H](/C(=C/C)/C)[C@@H](C)C[C@H]3O)O[C@@H](C2)CC=C(C[C@H](C=CC=C4[C@]5([C@H]1C=C([C@@H](O)C5OC4)COC(=O)CC6=COC=C6)O)C)C '
               'NAME: VM48642 REASON: CORRECT Found a spiroketal with a carbon '
               'connecting two distinct ring oxygens\n'
               ' * SMILES: '
               'O1[C@]2(O[C@@H]([C@H](CC2)C)C)C[C@]3(OC(=O)[C@]4([C@@]5(O)[C@](OCC5=CC=C[C@@H](CC(=CC[C@@]1(C3)[H])C)C)(\\C(=N\\O)\\C(=C4)C)[H])[H])[H] '
               'NAME: Milbemycin A3 oxime REASON: CORRECT Found a spiroketal '
               'with a carbon connecting two distinct ring oxygens\n'
               ' * SMILES: '
               'CC[C@@H](O)C[C@H]1CCC[C@@]2(C[C@@H]3OC(=O)\\C=C\\[C@@](C)(O)[C@@H](O)[C@H](C)[C@@H](O)[C@H](O[C@H]4CC[C@@H]([C@H](C)O4)N(C)C)[C@@H](O)[C@](C)(O)CCCCC\\C=C\\[C@@H]4CC(C)(C)O[C@@]4(O)C[C@H](O2)[C@H]3C)O1 '
               'NAME: ossamycin REASON: CORRECT Found a spiroketal with a '
               'carbon connecting two distinct ring oxygens\n'
               ' * SMILES: '
               'CCC(C)[C@H]1O[C@@]2(C[C@@H]3C[C@@H](C\\C=C(C)\\[C@@H](O)[C@@H](C)\\C=C\\C=C4/CO[C@@H]5[C@H](O)C(C)=C[C@@H](C(=O)O3)[C@]45O)O2)C[C@H](O)[C@@H]1C '
               'NAME: Avermectin B2a aglycone REASON: CORRECT Found a '
               'spiroketal with a carbon connecting two distinct ring oxygens\n'
               ' * SMILES: '
               'O=C1OC2C[C@]3(O[C@@H](C2)CC=C(C[C@H](C=CC=C(C4=C1C=C(C)C(=C4)O)C)C)C)O[C@@H]([C@@H](C)CC3)C '
               'NAME: Milbemycin beta3 REASON: CORRECT Found a spiroketal with '
               'a carbon connecting two distinct ring oxygens\n'
               ' * SMILES: '
               'NC(=O)CC[C@@H]1NC[C@@]2(OC[C@@H](O)[C@@H](O)[C@@H]2O)OC1=O '
               'NAME: alpha-chrysopine REASON: CORRECT Found a spiroketal with '
               'a carbon connecting two distinct ring oxygens\n'
               ' * SMILES: '
               'O=C1O[C@@H]2C[C@]3(O[C@@H](C2)CC=C(C[C@H](C=CC=C4[C@]5([C@H]1C=C([C@@H](OC)[C@H]5OC4)CO)O)C)C)O[C@@H]([C@@H](C)CC3)C '
               'NAME: Milbemycin alpha24 REASON: CORRECT Found a spiroketal '
               'with a carbon connecting two distinct ring oxygens\n'
               ' * SMILES: '
               'O=C1O[C@@H]2C[C@]3(O[C@@H](C2)CC=C(C[C@H](C=CC=C(C4=C1C=C(C)C(=C4)O)CO)C)C)O[C@@H]([C@@H](C)CC3)CC '
               'NAME: Milbemycin beta14 REASON: CORRECT Found a spiroketal '
               'with a carbon connecting two distinct ring oxygens\n'
               ' * SMILES: '
               'O1C2(OC3CC(OC(=O)C4C(O)(CC(OC)C(=C4)C)C(CO)=CC=CC(CC(=CC3)C)C)C2)CCC(C1CC)C '
               'NAME: Milbemycin beta2 REASON: CORRECT Found a spiroketal with '
               'a carbon connecting two distinct ring oxygens\n'
               'False positives: SMILES: '
               'O[C@@H]1Cc2c(O)cc3O[C@@]4(C[C@@H](c5c(O)cc(O)cc5O4)c3c2O[C@@H]1c1ccc(O)cc1)c1ccc(O)cc1 '
               'NAME: mahuannin D REASON: WRONGLY CLASSIFIED Found a '
               'spiroketal with a carbon connecting two distinct ring oxygens\n'
               ' * SMILES: O1C2(OCCCC2)CCC1 NAME: 1,6-Dioxaspiro[4.5]decane '
               'REASON: WRONGLY CLASSIFIED Found a spiroketal with a carbon '
               'connecting two distinct ring oxygens\n'
               ' * SMILES: '
               'O[C@@H]1Cc2c(O)cc(O)c([C@@H]3[C@@H](O)[C@H](Oc4c3c(O)cc3O[C@@]5(Oc6c([C@@H]7[C@@H](O)[C@H](Oc8cc(O)cc(O)c78)c7ccc(O)c(O)c7)c(O)cc(O)c6[C@@H]([C@H]5O)c43)c3ccc(O)c(O)c3)c3ccc(O)c(O)c3)c2O[C@@H]1c1ccc(O)c(O)c1 '
               'NAME: Cinnamtannin B2 REASON: WRONGLY CLASSIFIED Found a '
               'spiroketal with a carbon connecting two distinct ring oxygens\n'
               ' * SMILES: '
               'O1C2(OC(CC1(CC2)C)C3OC(CC3=O)(C(O)C4OC5(OC(CCC5)C(C)C(O)=O)CC4)C)C6OC(CC6)(CC(C)/C=C(/C=C/C7OC(C8(OCCC(C8O)C)O)CC7O)\\C)C '
               'NAME: Pectenotoxin 2 secoacid REASON: WRONGLY CLASSIFIED Found '
               'a spiroketal with a carbon connecting two distinct ring '
               'oxygens\n'
               ' * SMILES: '
               'O1C2=C(C=C(O)C3=C2[C@@]4(OC=5C=CC=C(C5[C@@]3(O)[C@H]4C)O)O[C@@H](C1(C)C)COC)C '
               'NAME: Cytorhizin C REASON: WRONGLY CLASSIFIED Found a '
               'spiroketal with a carbon connecting two distinct ring oxygens\n'
               ' * SMILES: '
               'C[C@@H]1C[C@H](OC(C)=O)[C@]2(COC(C)=O)[C@H](CC[C@H](O)[C@]22CO2)[C@@]1(C)[C@@H]1C[C@H]2CCO[C@H]2O1 '
               'NAME: ajubractin E REASON: WRONGLY CLASSIFIED Found a '
               'spiroketal with a carbon connecting two distinct ring oxygens\n'
               ' * SMILES: '
               'O=C1NC2=C(C(=O)C[C@H]3[C@]1([C@@]4([C@](O)(C5=CC(=O)[C@@H]6O[C@]5(OC6(C)C)CC4)CC3)C)C)C=C7CC8=C(C7=C2)CC(C)(C)OC8(C)C '
               'NAME: Shearinine H REASON: WRONGLY CLASSIFIED Found a '
               'spiroketal with a carbon connecting two distinct ring oxygens\n'
               ' * SMILES: '
               'O1[C@@]2([C@@]([C@@]3(C(C4C([C@@]5([C@@](CC4)(C[C@@H](O)[C@H](O)C5)[H])C)=CC3=O)C2)C)([C@@H]([C@]16OC[C@H](CC6)C)C)[H])[H] '
               'NAME: 9(11)-Dehydromanogenin REASON: WRONGLY CLASSIFIED Found '
               'a spiroketal with a carbon connecting two distinct ring '
               'oxygens\n'
               ' * SMILES: '
               'O1[C@]23O[C@H](C[C@](O)([C@@]3([C@@]4([C@]([C@]5([C@@]([C@@]6([C@@](CC5)(C([C@@H](O[C@@H]7OC[C@H](O)[C@H](O[C@@H]8O[C@@H]([C@@H](O)[C@H](O)[C@H]8O[C@@H]9OC[C@@H](O)[C@H](O)[C@H]9O)CO[C@@H]%10O[C@@H]([C@@H](O)[C@H](O)[C@H]%10O)CO)[C@H]7O[C@@H]%11O[C@H]([C@H](O)[C@@H](O)[C@H]%11O)C)CC6)(C)C)[H])C)(CC4)[H])C)(C2)C1)[H])[H])C)C=C(C)C '
               'NAME: Jujuboside A REASON: WRONGLY CLASSIFIED Found a '
               'spiroketal with a carbon connecting two distinct ring oxygens\n'
               ' * SMILES: '
               'O=C(O)C[C@H](OC)\\C=C/C=C\\C=C/C=C/C=C\\[C@@H]([C@H](O)[C@@H]([C@H]1O[C@]2(O[C@@H]([C@H](C)[C@@H](C2)OC)[C@@H]([C@@H](O)[C@H](C/C(=C/CC)/C)C)C)[C@@H](O)C[C@@H]1C)C)C '
               'NAME: Spirangien M REASON: WRONGLY CLASSIFIED Found a '
               'spiroketal with a carbon connecting two distinct ring oxygens\n'
               ' * SMILES: '
               'CO[C@H]1\\C=C\\O[C@@]2(C)OC3=C(C2=O)C2=C(C(O)=C3C)C([O-])=C(NC(=O)\\C(C)=C/C=C/[C@H](C)[C@H](O)[C@@H](C)[C@@H](O)[C@@H](C)[C@H](OC(C)=O)[C@@H]1C)\\C(=C/N(O)N1CCN(C)CC1)C2=O '
               "NAME: 2'-hydroxyrifampicin(1-) REASON: WRONGLY CLASSIFIED "
               'Found a spiroketal with a carbon connecting two distinct ring '
               'oxygens\n'
               ' * SMILES: '
               '[C@@H]1([C@@H]([C@H]([C@@H]([C@H](O1)CO)O)O)O)O[C@@]23[C@]4([C@@](C2)([C@@]5(O[C@]4(O[C@]3(C5)C)[H])O)[H])COC(C6=CC=C(C=C6)O)=O '
               'NAME: oxypaeoniflorin REASON: WRONGLY CLASSIFIED Found a '
               'spiroketal with a carbon connecting two distinct ring oxygens\n'
               ' * SMILES: '
               'O1C2C(C3(C(C4C(CC3)C5(C(=CC4)CC(O)CC5)C)C2)C)C(C16OC(CC6)(CO)CO)C '
               'NAME: Taccagenin REASON: WRONGLY CLASSIFIED Found a spiroketal '
               'with a carbon connecting two distinct ring oxygens\n'
               ' * SMILES: '
               'O=C1O[C@@H]([C@H](O)C=2C1=C(O)C=C3O[C@@]4(C[C@]5([C@]6([C@]4(CC23)O[C@](O5)(C(C)C)C6)C)O)C)C '
               'NAME: Talaromyolide A REASON: WRONGLY CLASSIFIED Found a '
               'spiroketal with a carbon connecting two distinct ring oxygens\n'
               ' * SMILES: '
               '[H][C@]12[C@@H](C)CCCCCCC[C@@]34O[C@@H]5[C@]6([H])[C@]7([H])O[C@]7(CO)[C@@H](O)[C@](O)(C(=O)[C@H]1C)[C@]2([H])[C@@]6(O3)[C@H](C)C[C@@]5(O4)C(C)=C '
               'NAME: wikstroelide E REASON: WRONGLY CLASSIFIED Found a '
               'spiroketal with a carbon connecting two distinct ring oxygens\n'
               ' * SMILES: '
               'C1(OC([C@]2([C@]3(CC1)[C@]4([C@H](CC2)C)O[C@@]5([C@@](C4)(C)[C@H]([C@@]6([C@@H]([C@](O5)(OC6=O)C)O)[H])C)OC3)[H])(C)C)=O '
               'NAME: fumigatonoid C REASON: WRONGLY CLASSIFIED Found a '
               'spiroketal with a carbon connecting two distinct ring oxygens\n'
               ' * SMILES: '
               'ClC1=C(O)C(Cl)=C2C(C=3C=C4CC5(OC6=C(C(C4=C(C3C(C2=C1O)=O)O)O5)C(=CC(=C6)OC)C)C)(C)C '
               'NAME: Bischloroanthrabenzoxocinone REASON: WRONGLY CLASSIFIED '
               'Found a spiroketal with a carbon connecting two distinct ring '
               'oxygens\n'
               ' * SMILES: '
               '[H][C@@]12O[C@]1(CO)[C@@H](O)[C@]1(O)C(=O)C(C)=C[C@@]1([H])[C@]13O[C@]4(O[C@H]([C@]21[H])[C@@](O4)([C@H](OC(C)=O)[C@H]3C)C(C)=C)\\C=C\\C=C\\CCCCCCCCC '
               'NAME: 12beta-acetoxyhuratoxin REASON: WRONGLY CLASSIFIED Found '
               'a spiroketal with a carbon connecting two distinct ring '
               'oxygens\n'
               ' * SMILES: O=C1C2=C(C(O)=C(C)C(=C2)O)CO[C@]13OCCC3 NAME: '
               'Terreinol REASON: WRONGLY CLASSIFIED Found a spiroketal with a '
               'carbon connecting two distinct ring oxygens\n'
               ' * SMILES: '
               'O1C2C(C3(C(C4C(CC3)C5(C(=CC4)CC(O)CC5)C)C2)C)C(C16OCC(O)(CC6)CO)C '
               'NAME: Spirotaccagenin REASON: WRONGLY CLASSIFIED Found a '
               'spiroketal with a carbon connecting two distinct ring oxygens\n'
               ' * SMILES: '
               'CC[C@H](C)[C@H]1O[C@@]2(C[C@@H]3C[C@@H](C\\C=C(C)\\[C@@H](O[C@H]4C[C@H](OC)[C@@H](O[C@H]5C[C@H](OC)[C@@H](O)[C@H](C)O5)[C@H](C)O4)[C@@H](C)\\C=C\\C=C4/CO[C@@H]5[C@H](O)C(C)=C[C@@H](C(=O)O3)[C@]45O)O2)C=C[C@@H]1C '
               'NAME: avermectin B1a REASON: WRONGLY CLASSIFIED Found a '
               'spiroketal with a carbon connecting two distinct ring oxygens\n'
               ' * SMILES: '
               'O1C2C(OC3(OCCC3)C(C2(O)C)C)C(O)C(C4OC5C(OC6C(OC7(C(OC8C(OC9C(OC%10C(OC%11CC%12OC%13C(OC(C=CC%13)/C=C\\C(O)CO)C(O)C%12OC%11C=C%10)C=C9)CC=CC8)CC7O)C6)C)CC(C5)C)CC14)C '
               'NAME: 54-Deoxy-50-hydroxyciguatoxin REASON: WRONGLY CLASSIFIED '
               'Found a spiroketal with a carbon connecting two distinct ring '
               'oxygens\n'
               ' * SMILES: '
               'O1[C@@]2([C@@]([C@@]3([C@]([C@]4([C@@]([C@@]5([C@@](CC4)(C[C@@H](OC(=O)C)CC5)[H])C)(CC3=O)[H])[H])(C2)[H])C)([C@@H]([C@]16OC[C@@H](CC6)C)C)[H])[H] '
               'NAME: HECOGENIN ACETATE REASON: WRONGLY CLASSIFIED Found a '
               'spiroketal with a carbon connecting two distinct ring oxygens\n'
               ' * SMILES: '
               'C[C@H]1[C@H]2[C@H](C[C@H]3[C@@H]4CC[C@@H]5C[C@@H](O)C[C@@H](O)[C@]5(C)[C@H]4CC[C@]23C)O[C@]11CCC(=C)CO1 '
               'NAME: Convallamarogenin REASON: WRONGLY CLASSIFIED Found a '
               'spiroketal with a carbon connecting two distinct ring oxygens\n'
               ' * SMILES: O1C(=CC(C(O)C)=C2C1OC(C)=C2OC)C NAME: Chenopodolan '
               'B REASON: WRONGLY CLASSIFIED Found a spiroketal with a carbon '
               'connecting two distinct ring oxygens\n'
               'False negatives: SMILES: '
               '[H][C@@]12CCC3=CC(=O)C=C[C@]3(C)[C@@]1(F)[C@@H](O)C[C@@]1(C)[C@@]2([H])C[C@H]2OC3(CCCC3)O[C@@]12C(=O)COC(C)=O '
               'NAME: amcinonide REASON: MISSED No spiroketal structure found\n'
               ' * SMILES: NC(=N)NCC1COC2(CCCCC2)O1 NAME: guanadrel REASON: '
               'MISSED No spiroketal structure found\n'
               ' * SMILES: CCCN(CC)CC1COC2(CCC(CC2)C(C)(C)C)O1 NAME: '
               'spiroxamine REASON: MISSED No spiroketal structure found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'O([C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)CO)[C@H]3[C@@H](O)[C@H](OC(O)[C@@H]3O)CO',
                                     'name': 'beta-D-Glcp-(1->2)-beta-D-Glcp-(1->3)-D-Galp',
                                     'reason': 'No spiroketal structure found'},
                                 {   'smiles': 'CCS(=O)(=O)NCC[C@@H]1CC[C@@H]([C@H](O1)CO)NC(=O)NC2=CC(=CC(=C2)Cl)Cl',
                                     'name': '1-(3,5-dichlorophenyl)-3-[(2S,3S,6S)-6-[2-(ethylsulfonylamino)ethyl]-2-(hydroxymethyl)-3-oxanyl]urea',
                                     'reason': 'No spiroketal structure found'},
                                 {   'smiles': 'C(C(O)=O)C/C=C\\C/C=C\\C\\C=C/C=C/C=C/[C@H]1[C@H](C/C=C\\CC)O1',
                                     'name': '(16S,17S)-epoxy-(4Z,7Z,10Z,12E,14E,19Z)-docosahexaenoic '
                                             'acid',
                                     'reason': 'No spiroketal structure found'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@H]2[C@@H](CO)O[C@@H](O[C@@H]3[C@@H](CO)O[C@@H](O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O',
                                     'name': 'alpha-D-Galp-(1->4)-beta-D-Galp-(1->4)-beta-D-Glcp',
                                     'reason': 'No spiroketal structure found'},
                                 {   'smiles': '[Li+].[Br-]',
                                     'name': 'lithium bromide',
                                     'reason': 'No spiroketal structure found'},
                                 {   'smiles': 'C=1(OC)C2=C(C=C3C[C@H]([C@](CC=4C=C(OC)C(OC)=C(C4C13)OC)(C)O)C)OCO2',
                                     'name': 'Besigomsin',
                                     'reason': 'No spiroketal structure found'},
                                 {   'smiles': 'C1=CC=C(C(=C1)C=CC(=O)C2=CC=CN2)Cl',
                                     'name': '3-(2-chlorophenyl)-1-(1H-pyrrol-2-yl)-2-propen-1-one',
                                     'reason': 'No spiroketal structure found'},
                                 {   'smiles': 'O=C1OC(O)C(=C1C(O)C(C)C)C',
                                     'name': '5-Hydroxy-3-(1-hydroxy-2-methylpropyl)-4-methyl-2(5H)-furanone',
                                     'reason': 'No spiroketal structure found'},
                                 {   'smiles': '[H][C@@]1(COC(C)(C)[C@@]1([H])OC(=O)\\C=C\\c1ccccc1)c1c(O)ccc2C(=O)C[C@H](Oc12)c1ccccc1',
                                     'name': '(+)-tephrorin B',
                                     'reason': 'No spiroketal structure found'},
                                 {   'smiles': 'C[C@@H]1CN(C(=O)C2=C(C3=CC=CC=C3CO[C@@H]1CN(C)CC4=CC=CC(=C4)C(=O)O)C5=CC=CC=C5N2C)[C@@H](C)CO',
                                     'name': 'LSM-9341',
                                     'reason': 'No spiroketal structure '
                                               'found'}],
    'sample_false_negatives': [],
    'attempt': 5,
    'reasoning': 'To improve the previous attempt at classifying spiroketal '
                 "structures, let's analyze the previous code and its "
                 'outcomes. \n'
                 '\n'
                 'The key requirement for a spiroketal is that it contains a '
                 'spiro center, typically a carbon, connecting two distinct '
                 'rings each with an oxygen atom. The previous code tries to '
                 'detect this structure using a SMARTS pattern for a carbon '
                 'attached to two oxygens, but the pattern may miss some '
                 'structural nuances or ring associations. This might explain '
                 'some of the false positive and false negative '
                 'classifications.\n'
                 '\n'
                 'Here are some improvements and considerations for the '
                 'revised code:\n'
                 '\n'
                 '1. Ensure that the SMARTS pattern accurately captures the '
                 'definition of a spiroketal: a single carbon is shared by two '
                 'otherwise separate rings, and each ring contains one of the '
                 'attached oxygens.\n'
                 '\n'
                 '2. Verify that the found carbon is indeed part of a fused '
                 'bicyclic structure where it connects two oxygen-containing '
                 'rings, following the classic spiro configuration.\n'
                 '\n'
                 '3. Correct any logic-related issues that might result from '
                 'how substructure matching is conducted, ensuring the rings '
                 'are distinct and consistent with the definition.\n'
                 '\n'
                 '4. Enhance robustness in handling complex spiro structures '
                 'which might be missed by overly-simplified SMARTS '
                 'patterning.\n'
                 '\n'
                 "Here's how the revised program might look:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 101,
    'num_false_positives': 100,
    'num_true_negatives': 1425,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.5024875621890548,
    'recall': 1.0,
    'f1': 0.6688741721854305,
    'accuracy': 0.9384993849938499,
    'negative_predictive_value': 1.0}