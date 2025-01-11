"""
Classifies: CHEBI:28765 phosphatidylinositol phosphate
"""
"""
Classifies: phosphatidylinositol phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylinositol_phosphate(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol phosphate based on its SMILES string.
    PIPs contain a myo-inositol ring with at least one phosphate group attached directly
    to the ring hydroxyl groups (not counting the phosphate that links to glycerol).
    Most common forms are PI(3)P, PI(4)P, PI(5)P, PI(3,4)P2, PI(3,5)P2, PI(4,5)P2, and PI(3,4,5)P3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a PIP, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for inositol ring (6-membered ring with all carbons and correct substitution)
    inositol_pattern = Chem.MolFromSmarts("[C]1[C][C][C][C][C]1")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring found"

    # Find inositol ring atoms
    inositol_matches = mol.GetSubstructMatches(inositol_pattern)
    if not inositol_matches:
        return False, "No inositol ring found"
    inositol_atoms = set(inositol_matches[0])

    # Pattern for phosphate group attached to oxygen
    phosphate_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX2])[OX2]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate groups found"

    # Find all phosphates attached to inositol ring
    inositol_phosphates = 0
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    
    for match in phosphate_matches:
        o_atom = mol.GetAtomWithIdx(match[0])
        # Get the atom the oxygen is connected to
        for neighbor in o_atom.GetNeighbors():
            if neighbor.GetIdx() in inositol_atoms:
                inositol_phosphates += 1
                break

    # Must have at least one phosphate directly on inositol ring
    # (not counting the bridging phosphate to glycerol)
    if inositol_phosphates < 2:  # One phosphate is the bridge to glycerol
        return False, "No additional phosphate groups on inositol ring beyond the glycerol linkage"

    # Look for glycerol backbone pattern
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    
    # Verify either glycerol backbone or modified structure
    is_modified = False
    if not mol.HasSubstructMatch(glycerol_pattern):
        is_modified = True
        if not mol.HasSubstructMatch(ester_pattern):
            return False, "Neither glycerol backbone nor modified structure found"

    # For standard PIPs, verify fatty acid chains
    if not is_modified:
        fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
        fatty_acid_matches = len(mol.GetSubstructMatches(fatty_acid_pattern))
        if fatty_acid_matches < 2:
            return False, "Missing fatty acid chains in glycerol-based PIP"

    # Classify based on number of phosphates on inositol
    actual_modifications = inositol_phosphates - 1  # Subtract the bridging phosphate
    
    if actual_modifications >= 3:
        return True, "Contains inositol ring with trisphosphate modification (PIP3)"
    elif actual_modifications == 2:
        return True, "Contains inositol ring with bisphosphate modification (PIP2)"
    elif actual_modifications == 1:
        return True, "Contains inositol ring with monophosphate modification (PIP)"
    else:
        return False, "Insufficient phosphate modifications on inositol ring"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:28765',
                          'name': 'phosphatidylinositol phosphate',
                          'definition': 'Any member of the phosphoinositide '
                                        'family of compounds, of which seven '
                                        'occur naturally.',
                          'parents': ['CHEBI:18179'],
                          'xrefs': ['PMID:16441841', 'PMID:2538726'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.35064935064935066 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](OP(O)(O)=O)[C@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC '
               'NAME: '
               '1-stearoyl-2-arachidonoyl-sn-glycero-3-phospho-1D-myo-inositol '
               '4,5-biphosphate REASON: CORRECT Contains inositol ring with '
               'bisphosphate modification (PIP2)\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC '
               'NAME: '
               '1-stearoyl-2-arachidonoyl-sn-glycero-3-phospho-1D-myo-inositol '
               '5-phosphate REASON: CORRECT Contains inositol ring with '
               'monophosphate modification (PIP)\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](OP(O)(O)=O)[C@H](O)[C@H]1O)OC(=O)CCCCCCC\\C=C/CCCCCCCC '
               'NAME: 1-stearoyl-2-oleoyl-sn-glycero-3-phospho-1D-myo-inositol '
               '4-phosphate REASON: CORRECT Contains inositol ring with '
               'monophosphate modification (PIP)\n'
               ' * SMILES: '
               'CCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCC '
               'NAME: '
               '1,2-dioctanoyl-sn-glycero-3-phospho-(1D-myo-inositol-3,4-bisphosphate) '
               'REASON: CORRECT Contains inositol ring with bisphosphate '
               'modification (PIP2)\n'
               ' * SMILES: '
               '[C@@H]1([C@@H]([C@@H](O)[C@H]([C@H]([C@@H]1O)OP(=O)(OCOC(=O)CCC)OCOC(=O)CCC)O)OP(=O)(OCOC(=O)CCC)OCOC(=O)CCC)OP(=O)(OCOC(=O)CCC)OCOC(=O)CCC '
               'NAME: D-myo-Ins(1,4,5)P3 hexakis(butyryloxymethyl) ester '
               'REASON: CORRECT Contains inositol ring with bisphosphate '
               'modification (PIP2)\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCC\\C=C/CCCCCCCC '
               'NAME: 1-stearoyl-2-oleoyl-sn-glycero-3-phospho-1D-myo-inositol '
               '5-phosphate REASON: CORRECT Contains inositol ring with '
               'monophosphate modification (PIP)\n'
               ' * SMILES: '
               'CCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCC '
               'NAME: '
               "1,2-dioctanoyl-sn-glycero-3-phospho-(1'-myo-inositol-5'-phosphate) "
               'REASON: CORRECT Contains inositol ring with monophosphate '
               'modification (PIP)\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCCCCCCCCCC '
               'NAME: '
               '1,2-dihexadecanoyl-sn-glycero-3-phospho-(1D-myo-inositol-3,4,5-trisphosphate) '
               'REASON: CORRECT Contains inositol ring with trisphosphate '
               'modification (PIP3)\n'
               ' * SMILES: '
               'O=C(CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)O[C@@H](COP(O[C@H]1[C@@H]([C@H]([C@@H]([C@H]([C@H]1O)O)OP(O)(=O)O)O)O)(=O)O)COC(=O)CCCCCCCCCCCCCCCCC '
               'NAME: '
               '1-stearoyl-2-arachidonoyl-sn-glycero-3-phospho-1D-myo-inositol '
               '4-phosphate REASON: CORRECT Contains inositol ring with '
               'monophosphate modification (PIP)\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](OP(O)(O)=O)[C@H](O)[C@@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCCCCCCCCCC '
               'NAME: '
               "1,2-dipalmitoyl-sn-glycero-3-phospho-(1'D-myo-inositol-3',5'-bisphosphate) "
               'REASON: CORRECT Contains inositol ring with bisphosphate '
               'modification (PIP2)\n'
               ' * SMILES: '
               '[C@@H]1(C(C(C([C@H](C1O)OP(=O)(O)O)OP(=O)(O)O)O)O)OP(OC[C@@](COC(CCCCCCCCCCCCCCC)=O)(OC(CCCCCCCCCCCCCCCCC)=O)[H])(O)=O '
               "NAME: PIP2[3',4'](16:0/18:0) REASON: CORRECT Contains inositol "
               'ring with bisphosphate modification (PIP2)\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](OP(O)(O)=O)[C@H](O)[C@H]1O)OC(=O)CCCCCCCCCCCCCCC '
               'NAME: '
               '1,2-dihexadecanoyl-sn-glycero-3-phospho-(1D-myo-inositol-4-phosphate) '
               'REASON: CORRECT Contains inositol ring with monophosphate '
               'modification (PIP)\n'
               ' * SMILES: '
               'CCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCC '
               'NAME: '
               "1,2-dioctanoyl-sn-glycero-3-phospho-(1D-myo-inositol-3'-phosphate) "
               'REASON: CORRECT Contains inositol ring with monophosphate '
               'modification (PIP)\n'
               ' * SMILES: '
               '[C@@H]1([C@@H]([C@@H]([C@@H]([C@H]([C@@H]1O)OP(O)(=O)O)O)O)OP(OC[C@@H](COC(=O)CCCCCCC)OC(=O)CCCCCCC)(=O)O)O '
               'NAME: '
               '1,2-dioctanoyl-sn-glycero-3-phospho-(1D-myo-inositol-4-phosphate) '
               'REASON: CORRECT Contains inositol ring with monophosphate '
               'modification (PIP)\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCCCCCCCCCC '
               'NAME: '
               '1,2-dipalmitoyl-sn-glycero-3-phospho-(1D-myo-inositol-5-phosphate) '
               'REASON: CORRECT Contains inositol ring with monophosphate '
               'modification (PIP)\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCCCCCCCCCC '
               'NAME: '
               "1,2-dipalmitoyl-sn-glycero-3-phospho-(1'D-myo-inositol-3'-phosphate) "
               'REASON: CORRECT Contains inositol ring with monophosphate '
               'modification (PIP)\n'
               ' * SMILES: '
               '[C@@H]1(C(C(C([C@H](C1O)OP(=O)(O)O)O)O)O)OP(OC[C@](COC(CCCCCCCCCCCCCCCCC)=O)([H])OC(CCCCCCCCCCCCCCC)=O)(O)=O '
               'NAME: PIP(18:0/16:0) REASON: CORRECT Contains inositol ring '
               'with monophosphate modification (PIP)\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](OP(O)(O)=O)[C@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCCCCCCCCCC '
               'NAME: '
               '1,2-dihexadecanoyl-sn-glycero-3-phospho-(1D-myo-inositol-4,5-bisphosphate) '
               'REASON: CORRECT Contains inositol ring with bisphosphate '
               'modification (PIP2)\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCCCCCCCCCC '
               'NAME: '
               '1,2-dihexadecanoyl-sn-glycero-3-phospho-(1D-myo-inositol-3,4-bisphosphate) '
               'REASON: CORRECT Contains inositol ring with bisphosphate '
               'modification (PIP2)\n'
               ' * SMILES: '
               'CCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](OP(O)(O)=O)[C@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCC '
               'NAME: '
               '1,2-dioctanoyl-sn-glycero-3-phospho-(1D-myo-inositol-4,5-bisphosphate) '
               'REASON: CORRECT Contains inositol ring with bisphosphate '
               'modification (PIP2)\n'
               ' * SMILES: '
               'CCCCC\\C=C/C\\C=C/CCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCC\\C=C/C\\C=C/CCCCC '
               'NAME: 1,2-dilinoleoyl-sn-glycero-3-phospho-1D-myo-inositol '
               '5-phosphate REASON: CORRECT Contains inositol ring with '
               'monophosphate modification (PIP)\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC '
               'NAME: '
               '1-stearoyl-2-arachidonoyl-sn-glycero-3-phospho-(1D-myo-inositol '
               '3,4,5-triphosphate) REASON: CORRECT Contains inositol ring '
               'with trisphosphate modification (PIP3)\n'
               ' * SMILES: '
               'CCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCC '
               'NAME: '
               '1,2-dioctanoyl-sn-glycero-3-phospho-(1D-myo-inositol-3,4,5-trisphosphate) '
               'REASON: CORRECT Contains inositol ring with trisphosphate '
               'modification (PIP3)\n'
               ' * SMILES: '
               'O(C[C@](COP(O)(=O)O[C@H]1C(C(C([C@H](C1O)OP(=O)(O)O)O)O)O)(OC(CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)=O)[H])C(CCCCCCCCC/C=C\\CCCCCC)=O '
               'NAME: PIP(18:1(11Z)/20:4(5Z,8Z,11Z,14Z)) REASON: CORRECT '
               'Contains inositol ring with monophosphate modification (PIP)\n'
               ' * SMILES: '
               'CCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](OP(O)(O)=O)[C@H](O)[C@@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCC '
               'NAME: '
               '1,2-dioctanoyl-sn-glycero-3-phospho-(1D-myo-inositol-3,5-bisphosphate) '
               'REASON: CORRECT Contains inositol ring with bisphosphate '
               'modification (PIP2)\n'
               'False positives: SMILES: '
               'P(OC1C(O)C(O)C(O)[C@@H](O)C1O)(OC[C@H](OC(=O)CCCCCCC/C=C\\CCCCCCC)COC(=O)CCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)(O)=O '
               'NAME: PI(22:4(7Z,10Z,13Z,16Z)/17:1(9Z)) REASON: WRONGLY '
               'CLASSIFIED Contains inositol ring with monophosphate '
               'modification (PIP)\n'
               ' * SMILES: '
               'P(OC1C(O)C(O)C(O)[C@@H](O)C1O)(OC[C@H](OC(=O)CCCCCCCC(O)/C=C/C(O)=O)COC(=O)CCCCCCC/C=C\\CCCCCCCC)(O)=O '
               'NAME: OHDdiA-PI REASON: WRONGLY CLASSIFIED Contains inositol '
               'ring with monophosphate modification (PIP)\n'
               ' * SMILES: '
               'P(OC1C(O)C(O)C(O)[C@@H](O)C1O)(OC[C@H](OC(=O)CCCCCCCCCCCCCCCCCCCC)COC(=O)CCCCCCCCCCCCCCCC)(O)=O '
               'NAME: PI(17:0/21:0) REASON: WRONGLY CLASSIFIED Contains '
               'inositol ring with monophosphate modification (PIP)\n'
               ' * SMILES: '
               '[C@@H]1(C(C(C([C@H](C1O)O)O)O)O)OP(OC[C@@](COC(CCCCCCCCC/C=C\\CCCCCC)=O)(OC(CCCCCCCCC/C=C\\CCCCCCCC)=O)[H])(O)=O '
               'NAME: PI(18:1(11Z)/20:1(11Z)) REASON: WRONGLY CLASSIFIED '
               'Contains inositol ring with monophosphate modification (PIP)\n'
               ' * SMILES: '
               '[C@@H]1(C(C(C([C@H](C1O)O)O)O)O)OP(OC[C@@](COC(CCCCCC/C=C\\C/C=C\\C/C=C\\CCCCC)=O)(OC(CCCCCCCCC/C=C\\CCCCCC)=O)[H])(O)=O '
               'NAME: PI(20:3(8Z,11Z,14Z)/18:1(11Z)) REASON: WRONGLY '
               'CLASSIFIED Contains inositol ring with monophosphate '
               'modification (PIP)\n'
               ' * SMILES: '
               'P(OC1C(O)C(O)C(O)[C@@H](O)C1O)(OC[C@H](OC(=O)CCCCCCCCCCCC)COC(=O)CCCC/C=C\\C/C=C\\C/C=C\\CCCCC)(O)=O '
               'NAME: PI(18:3(6Z,9Z,12Z)/13:0) REASON: WRONGLY CLASSIFIED '
               'Contains inositol ring with monophosphate modification (PIP)\n'
               ' * SMILES: '
               'P(OC1C(O)C(O)C(O)[C@@H](O)C1O)(OC[C@H](OC(=O)CCCCCCC/C=C\\CCCCCC)CO/C=C\\CCCCCCCCCCCCCC)(O)=O '
               'NAME: PI(P-16:0/16:1(9Z)) REASON: WRONGLY CLASSIFIED Contains '
               'inositol ring with monophosphate modification (PIP)\n'
               ' * SMILES: '
               'P(OC1C(O)C(O)C(O)[C@@H](O)C1O)(OC[C@H](O)COC(=O)CCCCCCC/C=C\\CCCC)(O)=O '
               'NAME: PI(14:1(9Z)/0:0) REASON: WRONGLY CLASSIFIED Contains '
               'inositol ring with monophosphate modification (PIP)\n'
               ' * SMILES: '
               'P(OC1C(O)C(O)C(O)[C@@H](O)C1O)(OC[C@H](OC(=O)CCCCCCCCCCCCCCCCCCC)COC(=O)CCCCCCCCC/C=C\\CCCCCCCCCC)(O)=O '
               'NAME: PI(22:1(11Z)/20:0) REASON: WRONGLY CLASSIFIED Contains '
               'inositol ring with monophosphate modification (PIP)\n'
               ' * SMILES: '
               'P(OC1C(O)C(O)C(O)[C@@H](O)C1O)(OC[C@H](OC(=O)CCCCCCC/C=C\\CCCCCC)COC(=O)CCCCCCCCCCCCCC)(O)=O '
               'NAME: PI(15:0/16:1(9Z)) REASON: WRONGLY CLASSIFIED Contains '
               'inositol ring with monophosphate modification (PIP)\n'
               ' * SMILES: '
               'P(O[C@@H]1C(O)[C@H](OP(O)(O)=O)C(O)C(O)C1O)(OC[C@H](OC(=O)CCCCCCC/C=C\\CCCCCC)COC(=O)CCCCCCC/C=C\\CCCCCC)(O)=O '
               'NAME: PIP(16:1(9Z)/16:1(9Z)) REASON: WRONGLY CLASSIFIED '
               'Contains inositol ring with monophosphate modification (PIP)\n'
               ' * SMILES: '
               'P(OC1C(O)C(O)C(O)[C@@H](O)C1O)(OC[C@H](OC(=O)CCCCCCCCCCCCCCCC)COC(=O)CCCCCCCCCCCCCCCCCCCC)(O)=O '
               'NAME: PI(21:0/17:0) REASON: WRONGLY CLASSIFIED Contains '
               'inositol ring with monophosphate modification (PIP)\n'
               ' * SMILES: '
               'P(OC1C(O)C(O)C(O)[C@@H](O)C1O)(OC[C@H](OC(=O)CCCCCCCCCCCCCCCCCCCC)COC(=O)CCCCCCCCCCCCCCC)(O)=O '
               'NAME: PI(16:0/21:0) REASON: WRONGLY CLASSIFIED Contains '
               'inositol ring with monophosphate modification (PIP)\n'
               ' * SMILES: '
               '[C@]1(C(C(C([C@](C1(O[H])[H])(O[H])[H])(O[H])[H])(O[H])[H])(O[H])[H])(OP(OC([C@@](C(OC(C(C(C(C(C(C(C(C(C(/C(=C(\\C(C(C(C(C(C(C(C([H])([H])[H])([H])[H])([H])[H])([H])[H])([H])[H])([H])[H])([H])[H])([H])[H])/[H])/[H])([H])[H])([H])[H])([H])[H])([H])[H])([H])[H])([H])[H])([H])[H])([H])[H])([H])[H])=O)([H])[H])(OC(C(C(C(C(C(C(/C(=C(\\C(/C(=C(\\C(/C(=C(\\C(/C(=C(\\C(C([H])([H])[H])([H])[H])/[H])/[H])([H])[H])/[H])/[H])([H])[H])/[H])/[H])([H])[H])/[H])/[H])([H])[H])([H])[H])([H])[H])([H])[H])([H])[H])([H])[H])=O)[H])([H])[H])(O[H])=O)[H] '
               'NAME: PI(20:1(11Z)/20:4(8Z,11Z,14Z,17Z)) REASON: WRONGLY '
               'CLASSIFIED Contains inositol ring with monophosphate '
               'modification (PIP)\n'
               ' * SMILES: '
               'P(OC1C(O)C(O)C(O)[C@@H](O)C1O)(OC[C@H](OC(=O)CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCCCCCCCC/C=C\\CCCCCCCCCC)(O)=O '
               'NAME: PI(22:1(11Z)/18:3(9Z,12Z,15Z)) REASON: WRONGLY '
               'CLASSIFIED Contains inositol ring with monophosphate '
               'modification (PIP)\n'
               ' * SMILES: '
               'P(OC1C(O)C(O)C(O)[C@@H](O)C1O)(OC[C@H](OC(=O)CCCCCCC/C=C\\C/C=C\\CCCC)COCCCCCCCCCCCCCCCCCCCC)(O)=O '
               'NAME: PI(O-20:0/17:2(9Z,12Z)) REASON: WRONGLY CLASSIFIED '
               'Contains inositol ring with monophosphate modification (PIP)\n'
               ' * SMILES: '
               'P(OC1C(O)C(O)C(O)[C@@H](O)C1O)(OC[C@H](OC(=O)CCCCCCC/C=C\\C/C=C\\CCCCC)COC(=O)CCCCCCCCCCCC)(O)=O '
               'NAME: PI(13:0/18:2(9Z,12Z)) REASON: WRONGLY CLASSIFIED '
               'Contains inositol ring with monophosphate modification (PIP)\n'
               ' * SMILES: '
               '[C@@H]1(C(C(C([C@H](C1O)O)O)O)O)OP(OC[C@@](COC(CCCCCCCCC/C=C\\CCCCCC)=O)(OC(CCCC/C=C\\C/C=C\\C/C=C\\CCCCC)=O)[H])(O)=O '
               'NAME: PI(18:1(11Z)/18:3(6Z,9Z,12Z)) REASON: WRONGLY CLASSIFIED '
               'Contains inositol ring with monophosphate modification (PIP)\n'
               ' * SMILES: '
               'P(OC1C(O[C@@H]2OC([C@H](O)C(O)[C@@H]2O)CO)C(O)C(O)C(O)[C@H]1O[C@H]3OC([C@@H](O)C(O)[C@H]3O)CO)(OC[C@H](O)COC(=O)CCCCCCC/C=C\\C/C=C\\CCCCC)(O)=O '
               'NAME: LPIM2(18:2(9Z,12Z)/0:0) REASON: WRONGLY CLASSIFIED '
               'Contains inositol ring with monophosphate modification (PIP)\n'
               ' * SMILES: '
               'P(OC1[C@H](O)[C@H](O)C(O)[C@H](O)[C@H]1O)(OC[C@H](OC(=O)CCCC(O)/C=C/C=C/C\\C=C\\CCCCCCCC)COC(=O)CCCCCCCCCCC/C=C\\C/C=C\\CCCCC)(O)=O '
               'NAME: PI(22:2(13Z,16Z)/20:3(6,8,11)-OH(5)) REASON: WRONGLY '
               'CLASSIFIED Contains inositol ring with monophosphate '
               'modification (PIP)\n'
               ' * SMILES: '
               '[C@@H]1(C(C(C([C@H](C1O)O)O)O)O)OP(OC[C@@](COC(CCCCCCC/C=C\\C/C=C\\CCCCC)=O)(OC(CCCCCCCCCCCCCCCCCCC)=O)[H])(O)=O '
               'NAME: PI(18:2(9Z,12Z)/20:0) REASON: WRONGLY CLASSIFIED '
               'Contains inositol ring with monophosphate modification (PIP)\n'
               ' * SMILES: '
               'P(OC1C(O[C@@H]2OC([C@H](O)C(O)[C@@H]2O)CO)C(O)C(O)C(O)[C@H]1O[C@H]3OC([C@@H](O)C(O)[C@H]3O)CO)(OC[C@H](O)COC(=O)CCCCCCCCCCCCCCCC)(O)=O '
               'NAME: LPIM2(17:0/0:0) REASON: WRONGLY CLASSIFIED Contains '
               'inositol ring with monophosphate modification (PIP)\n'
               ' * SMILES: '
               'P(OC1C(O)C(O)C(O)[C@@H](O)C1O)(OC[C@H](OC(=O)CCCCCCCC(=O)/C=C/C=O)COC(=O)CCCCCCC/C=C\\CCCCCCCC)(O)=O '
               'NAME: OKODA-PI REASON: WRONGLY CLASSIFIED Contains inositol '
               'ring with monophosphate modification (PIP)\n'
               ' * SMILES: '
               '[C@@H]1(C(C(C([C@H](C1O)O)O)O)O)OP(OC[C@@](COC(CCCCCCCCC/C=C\\C/C=C\\CCCCC)=O)(OC(CCCCCCCCCCCCCCCCC)=O)[H])(O)=O '
               'NAME: PI(20:2(11Z,14Z)/18:0) REASON: WRONGLY CLASSIFIED '
               'Contains inositol ring with monophosphate modification (PIP)\n'
               ' * SMILES: '
               '[C@@H]1(C(C(C([C@H](C1O)O)O)O)O)OP(OC[C@@](COC(CCCC/C=C\\C/C=C\\C/C=C\\CCCCC)=O)(OC(CCCCCCCCCCCCCCCCC)=O)[H])(O)=O '
               'NAME: PI(18:3(6Z,9Z,12Z)/18:0) REASON: WRONGLY CLASSIFIED '
               'Contains inositol ring with monophosphate modification (PIP)\n'
               'False negatives: NONE\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'Cc1cc(C)c2Cc3c(C)cc(C)c(Cc4c(C)cc(C)c(Cc5c(C)cc(C)c(Cc1c2C)c5C)c4C)c3C',
                                     'name': '4,6,10,12,16,18,22,24,25,26,27,28-dodecamethylcalix[4]arene',
                                     'reason': 'No inositol ring found'},
                                 {   'smiles': 'O=C1C(=C2C=C3[C@]([C@@H](C(C)C)[C@@H]([C@H]3O)OC(=O)C)(C)CC[C@]2(C)CC1)COC(=O)C',
                                     'name': 'Dahliane E',
                                     'reason': 'No inositol ring found'},
                                 {   'smiles': 'O[C@@H]([C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)N[C@@H](CC(C)C)C(O)=O)C',
                                     'name': 'Glu-Thr-Leu',
                                     'reason': 'No inositol ring found'},
                                 {   'smiles': 'CCOc1ccc(NC(=O)C(C)O)cc1',
                                     'name': 'p-Lactophenetide',
                                     'reason': 'No inositol ring found'},
                                 {   'smiles': 'O=C1NCC=2C1=C3C(N([C@@H]4O[C@H]([C@@H](O)[C@H]([C@H]4OC)O)C)C5=C3C=CC=C5)=C6NC7=C(C26)C=CC=C7',
                                     'name': "3'-epi-5'-methoxy-K252d",
                                     'reason': 'No inositol ring found'},
                                 {   'smiles': 'O=C(OC1C(O)C(OC(C1O)CO)OCC2OC(OCC(OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)C(O)C(C2O)O)CCCCCCCCCCCCCCC',
                                     'name': '[(2S)-2-hexadecanoyloxy-3-[(2S,3R,4S,5S,6R)-6-[[(2S,3R,4S,5S,6R)-4-hexadecanoyloxy-3,5-dihydroxy-6-(hydroxymethyl)tetrahydropyran-2-yl]oxymethyl]-3,4,5-trihydroxy-tetrahydropyran-2-yl]oxy-propyl] '
                                             '(9E,12E,15E)-octadeca-9,12,15-trienoate',
                                     'reason': 'No inositol ring found'},
                                 {   'smiles': 'O=C1C2=C(O)C=C(OC)C=C2C(=O)C3=C1[C@@H]([C@@H](O)[C@]([C@@H]3O)(O)C)C',
                                     'name': 'Altersolanol G',
                                     'reason': 'No inositol ring found'},
                                 {   'smiles': '[H][C@]1(O[C@](O)(C[C@H](O)[C@H]1NC(=O)CO)C([O-])=O)[C@H](O)[C@H](O)CO',
                                     'name': 'N-glycoloyl-alpha-neuraminate',
                                     'reason': 'No inositol ring found'},
                                 {   'smiles': 'OC(C(O)C/C=C\\C/C=C\\C/C=C\\CC)C/C=C\\C/C=C\\CCC(O)=O',
                                     'name': '10,11-DiHDPE',
                                     'reason': 'No inositol ring found'},
                                 {   'smiles': '[Na+].[H][C@]12SCC(C)=C(N1C(=O)[C@H]2NC(=O)[C@H](N)c1ccccc1)C([O-])=O',
                                     'name': 'cephalexin sodium',
                                     'reason': 'No inositol ring found'}],
    'sample_false_negatives': [   {   'smiles': 'CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC',
                                      'name': '1-stearoyl-2-arachidonoyl-sn-glycero-3-phospho-1D-myo-inositol '
                                              '5-phosphate',
                                      'reason': 'No additional phosphate '
                                                'groups on inositol ring '
                                                'beyond the glycerol linkage'},
                                  {   'smiles': 'CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](OP(O)(O)=O)[C@H](O)[C@H]1O)OC(=O)CCCCCCC\\C=C/CCCCCCCC',
                                      'name': '1-stearoyl-2-oleoyl-sn-glycero-3-phospho-1D-myo-inositol '
                                              '4-phosphate',
                                      'reason': 'No additional phosphate '
                                                'groups on inositol ring '
                                                'beyond the glycerol linkage'},
                                  {   'smiles': 'CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCC\\C=C/CCCCCCCC',
                                      'name': '1-stearoyl-2-oleoyl-sn-glycero-3-phospho-1D-myo-inositol '
                                              '5-phosphate',
                                      'reason': 'No additional phosphate '
                                                'groups on inositol ring '
                                                'beyond the glycerol linkage'},
                                  {   'smiles': 'CCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCC',
                                      'name': "1,2-dioctanoyl-sn-glycero-3-phospho-(1'-myo-inositol-5'-phosphate)",
                                      'reason': 'No additional phosphate '
                                                'groups on inositol ring '
                                                'beyond the glycerol linkage'},
                                  {   'smiles': 'O=C(CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)O[C@@H](COP(O[C@H]1[C@@H]([C@H]([C@@H]([C@H]([C@H]1O)O)OP(O)(=O)O)O)O)(=O)O)COC(=O)CCCCCCCCCCCCCCCCC',
                                      'name': '1-stearoyl-2-arachidonoyl-sn-glycero-3-phospho-1D-myo-inositol '
                                              '4-phosphate',
                                      'reason': 'No additional phosphate '
                                                'groups on inositol ring '
                                                'beyond the glycerol linkage'},
                                  {   'smiles': 'CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](OP(O)(O)=O)[C@H](O)[C@H]1O)OC(=O)CCCCCCCCCCCCCCC',
                                      'name': '1,2-dihexadecanoyl-sn-glycero-3-phospho-(1D-myo-inositol-4-phosphate)',
                                      'reason': 'No additional phosphate '
                                                'groups on inositol ring '
                                                'beyond the glycerol linkage'},
                                  {   'smiles': 'CCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCC',
                                      'name': "1,2-dioctanoyl-sn-glycero-3-phospho-(1D-myo-inositol-3'-phosphate)",
                                      'reason': 'No additional phosphate '
                                                'groups on inositol ring '
                                                'beyond the glycerol linkage'},
                                  {   'smiles': 'CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCCCCCCCCCC',
                                      'name': '1,2-dipalmitoyl-sn-glycero-3-phospho-(1D-myo-inositol-5-phosphate)',
                                      'reason': 'No additional phosphate '
                                                'groups on inositol ring '
                                                'beyond the glycerol linkage'},
                                  {   'smiles': 'CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCCCCCCCCCC',
                                      'name': "1,2-dipalmitoyl-sn-glycero-3-phospho-(1'D-myo-inositol-3'-phosphate)",
                                      'reason': 'No additional phosphate '
                                                'groups on inositol ring '
                                                'beyond the glycerol linkage'},
                                  {   'smiles': 'CCCCC\\C=C/C\\C=C/CCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCC\\C=C/C\\C=C/CCCCC',
                                      'name': '1,2-dilinoleoyl-sn-glycero-3-phospho-1D-myo-inositol '
                                              '5-phosphate',
                                      'reason': 'No additional phosphate '
                                                'groups on inositol ring '
                                                'beyond the glycerol linkage'}],
    'attempt': 2,
    'reasoning': 'Let me analyze the issues with the previous code:\n'
                 '\n'
                 '1. False Positives Analysis:\n'
                 '- The code incorrectly classified many phosphatidylinositols '
                 '(PIs) as PIPs\n'
                 "- The main issue is that it's not properly distinguishing "
                 'between a regular PI (with just one phosphate connecting '
                 'glycerol to inositol) and PIPs (which have additional '
                 'phosphates directly on the inositol ring)\n'
                 "- The phosphate counting was too simplistic and didn't "
                 'verify the position of phosphates\n'
                 '\n'
                 '2. Key Improvements Needed:\n'
                 '- We need to specifically identify phosphates directly '
                 'attached to the inositol ring hydroxyl groups\n'
                 '- We should exclude the bridging phosphate that connects to '
                 'glycerol when counting phosphate modifications\n'
                 '- The inositol ring pattern should be more specific to '
                 'myo-inositol configuration\n'
                 '- We should better validate the positions of phosphate '
                 'groups\n'
                 '\n'
                 "Here's my improved approach:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 16,
    'num_false_positives': 15,
    'num_true_negatives': 142258,
    'num_false_negatives': 11,
    'num_negatives': None,
    'precision': 0.5161290322580645,
    'recall': 0.5925925925925926,
    'f1': 0.5517241379310345,
    'accuracy': 0.9998172874209417,
    'negative_predictive_value': 0.9999226816804785}