"""
Classifies: CHEBI:18179 phosphoinositide
"""
"""
Classifies: phosphoinositide
Definition: Any phosphatidylinositol that is phosphorylated at one or more of the hydroxy groups of inositol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phosphoinositide(smiles: str):
    """
    Determines if a molecule is a phosphoinositide based on its SMILES string.
    Must have an inositol ring with at least one additional phosphate group beyond
    the one connecting to glycerol/lipid.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a phosphoinositide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for inositol ring core (cyclohexane with 6 oxygens)
    inositol_pattern = Chem.MolFromSmarts("[C]1([O])[C]([O])[C]([O])[C]([O])[C]([O])[C]1[O]")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring found"

    # Find all phosphate groups attached to the inositol ring
    # This pattern matches phosphates directly on the ring
    inositol_phosphate = Chem.MolFromSmarts("[C]1[C][C][C][C][C]1OP(=O)([O,OH])[O,OH]")
    phosphate_matches = mol.GetSubstructMatches(inositol_phosphate)
    
    if not phosphate_matches:
        return False, "No phosphate groups attached to inositol"

    # Count phosphates that are NOT connecting to a carbon chain (i.e., not the glycerol-connecting phosphate)
    additional_phosphate = Chem.MolFromSmarts("[C]1[C][C][C][C][C]1OP(=O)([O,OH])[O,OH]")
    phosphate_count = len(mol.GetSubstructMatches(additional_phosphate))
    
    # Look for glycerol/lipid attachment via phosphate
    glycerol_phosphate = Chem.MolFromSmarts("[CH2X4]-[CHX4]-[CH2X4]-OP(=O)([O,OH])-O")
    lipid_phosphate = Chem.MolFromSmarts("COP(=O)([O,OH])O[C]1[C][C][C][C][C]1")
    
    has_glycerol = mol.HasSubstructMatch(glycerol_phosphate)
    has_lipid = mol.HasSubstructMatch(lipid_phosphate)
    
    if not (has_glycerol or has_lipid):
        # Special case: Check for modified phosphoinositides without typical glycerol
        modified_phosphate_pattern = Chem.MolFromSmarts("[C]1[C][C][C][C][C]1OP(=O)(O)O")
        if mol.GetSubstructMatches(modified_phosphate_pattern):
            if phosphate_count >= 2:  # Need at least 2 phosphates for modified structures
                return True, f"Modified phosphoinositide with {phosphate_count} phosphate groups"
        return False, "No phospholipid attachment found"

    # For standard phosphoinositides, we need at least one additional phosphate
    # beyond the one connecting to glycerol/lipid
    if phosphate_count <= 1:
        return False, "Only has connecting phosphate, no additional phosphorylation"

    # Look for ester groups (fatty acid attachments) for standard phosphoinositides
    if has_glycerol:
        ester_pattern = Chem.MolFromSmarts("C(=O)OC")
        ester_matches = mol.GetSubstructMatches(ester_pattern)
        if len(ester_matches) < 1:
            return False, "No ester groups found"

    return True, f"Phosphoinositide with {phosphate_count} phosphate group(s) on inositol ring"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:18179',
                          'name': 'phosphoinositide',
                          'definition': 'Any phosphatidylinositol that is '
                                        'phosphorylated at one or more of the '
                                        'hydroxy groups of inositol.',
                          'parents': ['CHEBI:28874'],
                          'xrefs': [   'PMID:10782093',
                                       'PMID:15269334',
                                       'PMID:17035995',
                                       'PMID:19154715'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.3398692810457516 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               '[C@@H]1([C@@H]([C@@H]([C@@H]([C@H]([C@@H]1O)OP(O)(=O)O)O)O)OP(OC[C@@H](COC(=O)CCCCCCC)OC(=O)CCCCCCC)(=O)O)O '
               'NAME: '
               '1,2-dioctanoyl-sn-glycero-3-phospho-(1D-myo-inositol-4-phosphate) '
               'REASON: CORRECT Phosphoinositide with 2 phosphate group(s) on '
               'inositol ring\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC '
               'NAME: '
               '1-stearoyl-2-arachidonoyl-sn-glycero-3-phospho-1D-myo-inositol '
               '5-phosphate REASON: CORRECT Phosphoinositide with 2 phosphate '
               'group(s) on inositol ring\n'
               ' * SMILES: '
               '[C@@H]1(C(C(C([C@H](C1O)OP(=O)(O)O)OP(=O)(O)O)O)O)OP(OC[C@@](COC(CCCCCCCCCCCCCCC)=O)(OC(CCCCCCCCCCCCCCCCC)=O)[H])(O)=O '
               "NAME: PIP2[3',4'](16:0/18:0) REASON: CORRECT Phosphoinositide "
               'with 3 phosphate group(s) on inositol ring\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](OP(O)(O)=O)[C@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCCCCCCCCCC '
               'NAME: '
               '1,2-dihexadecanoyl-sn-glycero-3-phospho-(1D-myo-inositol-4,5-bisphosphate) '
               'REASON: CORRECT Phosphoinositide with 3 phosphate group(s) on '
               'inositol ring\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](OP(O)(O)=O)[C@H](O)[C@H]1O)OC(=O)CCCCCCC\\C=C/CCCCCCCC '
               'NAME: 1-stearoyl-2-oleoyl-sn-glycero-3-phospho-1D-myo-inositol '
               '4-phosphate REASON: CORRECT Phosphoinositide with 2 phosphate '
               'group(s) on inositol ring\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](OP(O)(O)=O)[C@H](O)[C@@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCCCCCCCCCC '
               'NAME: '
               "1,2-dipalmitoyl-sn-glycero-3-phospho-(1'D-myo-inositol-3',5'-bisphosphate) "
               'REASON: CORRECT Phosphoinositide with 3 phosphate group(s) on '
               'inositol ring\n'
               ' * SMILES: '
               '[C@@H]1(C(C(C([C@H](C1O)OP(=O)(O)O)O)O)O)OP(OC[C@](COC(CCCCCCCCCCCCCCCCC)=O)([H])OC(CCCCCCCCCCCCCCC)=O)(O)=O '
               'NAME: PIP(18:0/16:0) REASON: CORRECT Phosphoinositide with 2 '
               'phosphate group(s) on inositol ring\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](OP(O)(O)=O)[C@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC '
               'NAME: '
               '1-stearoyl-2-arachidonoyl-sn-glycero-3-phospho-1D-myo-inositol '
               '4,5-biphosphate REASON: CORRECT Phosphoinositide with 3 '
               'phosphate group(s) on inositol ring\n'
               ' * SMILES: '
               'CCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCC '
               'NAME: '
               '1,2-dioctanoyl-sn-glycero-3-phospho-(1D-myo-inositol-3,4,5-trisphosphate) '
               'REASON: CORRECT Phosphoinositide with 4 phosphate group(s) on '
               'inositol ring\n'
               ' * SMILES: '
               'P(O[C@@H]1C(O)[C@H](OP(O)(O)=O)C(OP(O)(O)=O)C(O)C1O)(OC[C@H](OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCCCCCCCCCC)(O)=O '
               'NAME: PIP2(16:0/16:0) REASON: CORRECT Phosphoinositide with 3 '
               'phosphate group(s) on inositol ring\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCCCCCCCCCC '
               'NAME: '
               '1,2-dipalmitoyl-sn-glycero-3-phospho-(1D-myo-inositol-5-phosphate) '
               'REASON: CORRECT Phosphoinositide with 2 phosphate group(s) on '
               'inositol ring\n'
               ' * SMILES: '
               'CCCCC\\C=C/C\\C=C/CCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCC\\C=C/C\\C=C/CCCCC '
               'NAME: 1,2-dilinoleoyl-sn-glycero-3-phospho-1D-myo-inositol '
               '5-phosphate REASON: CORRECT Phosphoinositide with 2 phosphate '
               'group(s) on inositol ring\n'
               ' * SMILES: '
               'CCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](OP(O)(O)=O)[C@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCC '
               'NAME: '
               '1,2-dioctanoyl-sn-glycero-3-phospho-(1D-myo-inositol-4,5-bisphosphate) '
               'REASON: CORRECT Phosphoinositide with 3 phosphate group(s) on '
               'inositol ring\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCCCCCCCCCC '
               'NAME: '
               "1,2-dipalmitoyl-sn-glycero-3-phospho-(1'D-myo-inositol-3'-phosphate) "
               'REASON: CORRECT Phosphoinositide with 2 phosphate group(s) on '
               'inositol ring\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](OP(O)(O)=O)[C@H](O)[C@H]1O)OC(=O)CCCCCCCCCCCCCCC '
               'NAME: '
               '1,2-dihexadecanoyl-sn-glycero-3-phospho-(1D-myo-inositol-4-phosphate) '
               'REASON: CORRECT Phosphoinositide with 2 phosphate group(s) on '
               'inositol ring\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC '
               'NAME: '
               '1-stearoyl-2-arachidonoyl-sn-glycero-3-phospho-(1D-myo-inositol '
               '3,4,5-triphosphate) REASON: CORRECT Phosphoinositide with 4 '
               'phosphate group(s) on inositol ring\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCC\\C=C/CCCCCCCC '
               'NAME: 1-stearoyl-2-oleoyl-sn-glycero-3-phospho-1D-myo-inositol '
               '5-phosphate REASON: CORRECT Phosphoinositide with 2 phosphate '
               'group(s) on inositol ring\n'
               ' * SMILES: '
               'CCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCC '
               'NAME: '
               '1,2-dioctanoyl-sn-glycero-3-phospho-(1D-myo-inositol-3,4-bisphosphate) '
               'REASON: CORRECT Phosphoinositide with 3 phosphate group(s) on '
               'inositol ring\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCCCCCCCCCC '
               'NAME: '
               '1,2-dihexadecanoyl-sn-glycero-3-phospho-(1D-myo-inositol-3,4,5-trisphosphate) '
               'REASON: CORRECT Phosphoinositide with 4 phosphate group(s) on '
               'inositol ring\n'
               ' * SMILES: '
               'CCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](OP(O)(O)=O)[C@H](O)[C@@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCC '
               'NAME: '
               '1,2-dioctanoyl-sn-glycero-3-phospho-(1D-myo-inositol-3,5-bisphosphate) '
               'REASON: CORRECT Phosphoinositide with 3 phosphate group(s) on '
               'inositol ring\n'
               ' * SMILES: '
               'CCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCC '
               'NAME: '
               "1,2-dioctanoyl-sn-glycero-3-phospho-(1D-myo-inositol-3'-phosphate) "
               'REASON: CORRECT Phosphoinositide with 2 phosphate group(s) on '
               'inositol ring\n'
               ' * SMILES: '
               'O=C(CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)O[C@@H](COP(O[C@H]1[C@@H]([C@H]([C@@H]([C@H]([C@H]1O)O)OP(O)(=O)O)O)O)(=O)O)COC(=O)CCCCCCCCCCCCCCCCC '
               'NAME: '
               '1-stearoyl-2-arachidonoyl-sn-glycero-3-phospho-1D-myo-inositol '
               '4-phosphate REASON: CORRECT Phosphoinositide with 2 phosphate '
               'group(s) on inositol ring\n'
               ' * SMILES: '
               'O(C[C@](COP(O)(=O)O[C@H]1C(C(C([C@H](C1O)OP(=O)(O)O)O)O)O)(OC(CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)=O)[H])C(CCCCCCCCC/C=C\\CCCCCC)=O '
               'NAME: PIP(18:1(11Z)/20:4(5Z,8Z,11Z,14Z)) REASON: CORRECT '
               'Phosphoinositide with 2 phosphate group(s) on inositol ring\n'
               ' * SMILES: '
               'P(O[C@@H]1C(O)[C@H](OP(O)(O)=O)C(OP(O)(O)=O)C(O)C1O)(OC[C@H](OC(=O)CCCCCCC/C=C\\C/C=C\\CCCCC)COC(=O)CCCCCCC/C=C\\C/C=C\\CCCCC)(O)=O '
               'NAME: PIP2(18:2/18:2) REASON: CORRECT Phosphoinositide with 3 '
               'phosphate group(s) on inositol ring\n'
               ' * SMILES: '
               'CCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCC '
               'NAME: '
               "1,2-dioctanoyl-sn-glycero-3-phospho-(1'-myo-inositol-5'-phosphate) "
               'REASON: CORRECT Phosphoinositide with 2 phosphate group(s) on '
               'inositol ring\n'
               'False positives: SMILES: '
               '[C@@H]1(C(C(C([C@H](C1O)O)O)O)O)OP(OC[C@@](COC(CCC/C=C\\C/C=C\\C/C=C\\CCCCCCCC)=O)(OC(CCC/C=C\\C/C=C\\C/C=C\\CCCCCCCC)=O)[H])(O)=O '
               'NAME: PI(20:3(5Z,8Z,11Z)/20:3(5Z,8Z,11Z)) REASON: WRONGLY '
               'CLASSIFIED Phosphoinositide with 1 phosphate group(s) on '
               'inositol ring\n'
               ' * SMILES: '
               'P(OC1C(O)C(O)C(O)[C@@H](O)C1O)(OC[C@H](O)COC(=O)CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)(O)=O '
               'NAME: PI(20:5(5Z,8Z,11Z,14Z,17Z)/0:0) REASON: WRONGLY '
               'CLASSIFIED Phosphoinositide with 1 phosphate group(s) on '
               'inositol ring\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP([O-])(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](OP([O-])([O-])=O)[C@H](O)[C@H]1O)OC(=O)CCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC '
               'NAME: '
               '1-stearoyl-2-arachidonoyl-sn-glycero-3-phospho-1D-myo-inositol '
               '4-phosphate(3-) REASON: WRONGLY CLASSIFIED Phosphoinositide '
               'with 2 phosphate group(s) on inositol ring\n'
               ' * SMILES: '
               'P(OC1C(O)C(O)C(O)[C@@H](O)C1O)(OC[C@H](OC(=O)CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCCCCCCCCCCCC)(O)=O '
               'NAME: PI(14:0/20:5(5Z,8Z,11Z,14Z,17Z)) REASON: WRONGLY '
               'CLASSIFIED Phosphoinositide with 1 phosphate group(s) on '
               'inositol ring\n'
               ' * SMILES: '
               '[C@@H]1(C(C(C([C@H](C1O)O)O)O)O)OP(OC[C@@](COC(CCCCCCCCCCCCCCC)=O)(OC(CCCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)=O)[H])(O)=O '
               'NAME: PI(16:0/20:4(8Z,11Z,14Z,17Z)) REASON: WRONGLY CLASSIFIED '
               'Phosphoinositide with 1 phosphate group(s) on inositol ring\n'
               ' * SMILES: '
               'CCCC(=O)OC[C@H](COP([O-])(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](OP([O-])([O-])=O)[C@H]1O)OC(=O)CCC '
               'NAME: '
               "1,2-dibutyryl-sn-glycero-3-phospho-(1'D-myo-inositol-5'-phosphate)(3-) "
               'REASON: WRONGLY CLASSIFIED Phosphoinositide with 2 phosphate '
               'group(s) on inositol ring\n'
               ' * SMILES: '
               '[C@@H]1(C(C(C([C@H](C1O)O)O)O)O)OP(OC[C@@](COC(CCCCCCC/C=C\\C/C=C\\CCCCC)=O)(OC(CCCCCCCCC/C=C\\CCCCCCCC)=O)[H])(O)=O '
               'NAME: PI(18:2(9Z,12Z)/20:1(11Z)) REASON: WRONGLY CLASSIFIED '
               'Phosphoinositide with 1 phosphate group(s) on inositol ring\n'
               ' * SMILES: '
               'P(OC1C(O)C(O)C(O)[C@@H](O)C1O)(OC[C@H](OC(=O)CCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)COC(=O)CCCCCCC/C=C\\C/C=C\\CCCCC)(O)=O '
               'NAME: PI(18:2(9Z,12Z)/22:4(7Z,10Z,13Z,16Z)) REASON: WRONGLY '
               'CLASSIFIED Phosphoinositide with 1 phosphate group(s) on '
               'inositol ring\n'
               ' * SMILES: '
               'P(OC1[C@H](O)[C@@H](O)C(O)[C@@H](O)[C@H]1O)(OC[C@H](OC(=O)CCC/C=C\\[C@@H](O)/C=C/C=C\\C/C=C\\C=C\\[C@@H](O)C/C=C\\CC)COC(=O)CCCCCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)(O)=O '
               'NAME: '
               'PI(22:4(10Z,13Z,16Z,19Z)/22:6(5Z,8E,10Z,13Z,15E,19Z)-2OH(7S, '
               '17S)) REASON: WRONGLY CLASSIFIED Phosphoinositide with 1 '
               'phosphate group(s) on inositol ring\n'
               ' * SMILES: '
               'P(OC1[C@H](O)[C@H](O)C(O)[C@H](O)[C@H]1O)(OC[C@H](OC(=O)CC/C=C\\C/C=C\\C\\C=C/C=C/C(O)C/C=C\\C/C=C\\CC)COC(=O)CCCCCCCCCCCCCCC)(O)=O '
               'NAME: PI(16:0/22:6(4Z,7Z,10Z,12E,16Z,19Z)-OH(14)) REASON: '
               'WRONGLY CLASSIFIED Phosphoinositide with 1 phosphate group(s) '
               'on inositol ring\n'
               ' * SMILES: '
               '[C@@H]1(C(C(C([C@H](C1O)O)O)O)O)OP(OC[C@@](COC(CCCCCCCCCCCCCCCCCCC)=O)(OC(CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)=O)[H])(O)=O '
               'NAME: PI(20:0/20:4(5Z,8Z,11Z,14Z)) REASON: WRONGLY CLASSIFIED '
               'Phosphoinositide with 1 phosphate group(s) on inositol ring\n'
               ' * SMILES: '
               '[C@@H]1(C(C(C([C@H](C1O)O)O)O)O)OP(OC[C@@](COC(CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)=O)(OC(CCCCCCCCC/C=C\\CCCCCC)=O)[H])(O)=O '
               'NAME: PI(20:4(5Z,8Z,11Z,14Z)/18:1(11Z)) REASON: WRONGLY '
               'CLASSIFIED Phosphoinositide with 1 phosphate group(s) on '
               'inositol ring\n'
               ' * SMILES: '
               'P(OC1C(O)C(O)C(O)[C@@H](O)C1O)(OC[C@H](O)COC(=O)CCCCCCCCCCCCCCCCCCCC)(O)=O '
               'NAME: PI(21:0/0:0) REASON: WRONGLY CLASSIFIED Phosphoinositide '
               'with 1 phosphate group(s) on inositol ring\n'
               ' * SMILES: '
               'P(OC1C(O)C(O)C(O)[C@@H](O)C1O)(OC[C@H](OC(=O)CCCCCCCCCCCCCCCCCCCC)COC(=O)CCCCCCCCCCCCCCCC)(O)=O '
               'NAME: PI(17:0/21:0) REASON: WRONGLY CLASSIFIED '
               'Phosphoinositide with 1 phosphate group(s) on inositol ring\n'
               ' * SMILES: '
               'P(OC1C(O)C(O)C(O)[C@@H](O)C1O)(OC[C@H](OC(=O)CCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)COC(=O)CCCCCCCCCCCCCCCCCC)(O)=O '
               'NAME: PI(19:0/22:4(7Z,10Z,13Z,16Z)) REASON: WRONGLY CLASSIFIED '
               'Phosphoinositide with 1 phosphate group(s) on inositol ring\n'
               ' * SMILES: '
               '[C@@H]1(C(C(C([C@H](C1O)O)O)O)O)OP(OC[C@@](COC(CCCCCCCCCCCCCCCCC)=O)(OC(CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)=O)[H])(O)=O '
               'NAME: PI(18:0/22:6(4Z,7Z,10Z,13Z,16Z,19Z)) REASON: WRONGLY '
               'CLASSIFIED Phosphoinositide with 1 phosphate group(s) on '
               'inositol ring\n'
               ' * SMILES: '
               'P(OC1[C@@H](O)[C@H](O)C(O)[C@H](O)[C@@H]1O)(OC[C@H](OC(=O)CCCCC)COC(=O)CCCCC)(O)=O '
               'NAME: PtdIns-(1,2-dihexanoyl) REASON: WRONGLY CLASSIFIED '
               'Phosphoinositide with 1 phosphate group(s) on inositol ring\n'
               ' * SMILES: '
               'P(OC1C(O)C(O)C(O)[C@@H](O)C1O)(OC[C@H](OC(=O)CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCCCCCCCCCCCCCCC)(O)=O '
               'NAME: PI(17:0/20:5(5Z,8Z,11Z,14Z,17Z)) REASON: WRONGLY '
               'CLASSIFIED Phosphoinositide with 1 phosphate group(s) on '
               'inositol ring\n'
               ' * SMILES: '
               'P(OC1C(O)C(O)C(O)C(O)C1O)(OCC(OC(=O)CCC/C=C\\CC2OC2C/C=C\\C/C=C\\CCCCC)COC(=O)CCCCCCCCCCCCCCCCC)(O)=O '
               'NAME: Pi(18:0/8,9-eet) REASON: WRONGLY CLASSIFIED '
               'Phosphoinositide with 1 phosphate group(s) on inositol ring\n'
               ' * SMILES: '
               'P(O[C@H]1[C@H](O[C@H]2OC([C@@H](O)C(O)[C@H]2O)CO)C(O)[C@H](O)C(O)C1O)(OC[C@H](OC(=O)CCCCCCCCCCCCC)COC(=O)CCCCCCCCCCCCCCC)(O)=O '
               'NAME: PIM1(16:0/14:0) REASON: WRONGLY CLASSIFIED '
               'Phosphoinositide with 1 phosphate group(s) on inositol ring\n'
               ' * SMILES: '
               '[C@@H]1(C(C(C([C@H](C1O)O)O)O)O)OP(OC[C@@](COC(CCCCCCCCC/C=C\\CCCCCC)=O)(OC(CCCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)=O)[H])(O)=O '
               'NAME: PI(18:1(11Z)/20:4(8Z,11Z,14Z,17Z)) REASON: WRONGLY '
               'CLASSIFIED Phosphoinositide with 1 phosphate group(s) on '
               'inositol ring\n'
               ' * SMILES: '
               '[C@@H]1(C(C(C([C@H](C1O)O)O)O)O)OP(OC[C@@](COC(CCCCCCC/C=C\\CCCCCCCC)=O)(OC(CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)=O)[H])(O)=O '
               'NAME: PI(18:1(9Z)/18:3(9Z,12Z,15Z)) REASON: WRONGLY CLASSIFIED '
               'Phosphoinositide with 1 phosphate group(s) on inositol ring\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCC(=O)OC[C@H](COP([O-])(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)OC(=O)CCCCCCCCCCCCCCC '
               'NAME: '
               '1,2-dihexadecanoyl-sn-glycero-3-phospho-D-myo-inositol(1-) '
               'REASON: WRONGLY CLASSIFIED Phosphoinositide with 1 phosphate '
               'group(s) on inositol ring\n'
               ' * SMILES: '
               '[C@@H]1(C(C(C([C@H](C1O)O)O)O)O)OP(OC[C@@](COC(CCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)=O)(OC(CCCCCCCCCCCCCCCCC)=O)[H])(O)=O '
               'NAME: PI(22:5(7Z,10Z,13Z,16Z,19Z)/18:0) REASON: WRONGLY '
               'CLASSIFIED Phosphoinositide with 1 phosphate group(s) on '
               'inositol ring\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP([O-])(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O '
               'NAME: 1-hexadecanoyl-sn--glycero-3-phospho-D-myo-inositol(1-) '
               'REASON: WRONGLY CLASSIFIED Phosphoinositide with 1 phosphate '
               'group(s) on inositol ring\n'
               'False negatives: SMILES: '
               '[C@@H]1([C@@H]([C@@H](O)[C@H]([C@H]([C@@H]1O)OP(=O)(OCOC(=O)CCC)OCOC(=O)CCC)O)OP(=O)(OCOC(=O)CCC)OCOC(=O)CCC)OP(=O)(OCOC(=O)CCC)OCOC(=O)CCC '
               'NAME: D-myo-Ins(1,4,5)P3 hexakis(butyryloxymethyl) ester '
               'REASON: MISSED No glycerol-phosphate backbone found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C[C@@H](C1=CC=CC=C1)NC(=O)C[C@H]2CC[C@@H]([C@@H](O2)CO)NC(=O)CN3CCOCC3',
                                     'name': '2-[(2R,5S,6R)-6-(hydroxymethyl)-5-[[2-(4-morpholinyl)-1-oxoethyl]amino]-2-oxanyl]-N-[(1S)-1-phenylethyl]acetamide',
                                     'reason': 'No inositol ring found'},
                                 {   'smiles': 'C[N+](C)(C)[C@@H](Cc1c[nH]c(n1)S(=O)C[C@H](NC(=O)CC[C@H]([NH3+])C([O-])=O)C([O-])=O)C([O-])=O',
                                     'name': 'N(alpha)-(L-gamma-glutamyl)-hercynyl-L-cysteine '
                                             'sulfoxide(1-)',
                                     'reason': 'No inositol ring found'},
                                 {   'smiles': 'O(C1=CC=2[C@]3([C@](N(CC3)C)(N(C2C=C1)C)[H])C)C(=O)N4CCC=5C(C4)=CC=CC5',
                                     'name': 'quilostigmine',
                                     'reason': 'No inositol ring found'},
                                 {   'smiles': 'O[C@@H]1[C@]23[C@@]4(N(C[C@@]([C@]2(C[C@@]4([C@]56[C@]3(CC(=O)[C@](C5)(C([C@H]6O)=C)[H])[H])[H])[H])(CC1)C)CC)[H]',
                                     'name': 'Bullatine G',
                                     'reason': 'No inositol ring found'},
                                 {   'smiles': 'NC1=NC=NC2=C1N=CN2[C@@H]3O[C@H](COP(=O)(O)O)[C@@H](OC(=O)[C@@H](N)CCC(O)=O)[C@H]3O',
                                     'name': "3'-L-glutamyl-AMP",
                                     'reason': 'No inositol ring found'},
                                 {   'smiles': 'O1C2(C(C3(C(C4(C(CC3OC(=O)C)C(OC(=O)CC4)(C)C)C)CC2)C)CC15C6N(C=7C5=CC=CC7)C(=O)C(N6)C)C',
                                     'name': 'Teraspiridole C_130091',
                                     'reason': 'No inositol ring found'},
                                 {   'smiles': 'O(C1[C@@H](OC(=O)C)C(O[C@@H](OC2=C(OC3=C(C2=O)C(O)=CC(O[C@@H]4OC([C@@H](O)[C@H](O)C4O)CO)=C3CC=C(C)C)C5=CC=C(OC)C=C5)[C@H]1O)C)[C@@H]6OC[C@@H](O)[C@H](OC(=O)C)C6O',
                                     'name': 'Sempervirenoside A',
                                     'reason': 'No inositol ring found'},
                                 {   'smiles': 'COC[C@]1(C(=O)C2CCN1CC2)CO',
                                     'name': '(2S)-2-(hydroxymethyl)-2-(methoxymethyl)-1-azabicyclo[2.2.2]octan-3-one',
                                     'reason': 'No inositol ring found'},
                                 {   'smiles': 'Oc1c(C2CC(Cc3ccccc23)c2ccc(OCc3ccc(cc3)C(F)(F)F)cc2)c(=O)oc2ccccc12',
                                     'name': 'Flocoumafen',
                                     'reason': 'No inositol ring found'},
                                 {   'smiles': 'O[C@H]1CC=2C(N(C=3C1=CC=CC3)C(=O)N)=CC=CC2',
                                     'name': '(S)-MHD',
                                     'reason': 'No inositol ring found'}],
    'sample_false_negatives': [],
    'attempt': 3,
    'reasoning': 'Let me analyze the issues with the previous code:\n'
                 '\n'
                 '1. False Positives Analysis:\n'
                 '- The code is incorrectly classifying simple '
                 'phosphatidylinositols (PIs) as phosphoinositides\n'
                 '- The main distinction is that phosphoinositides must have '
                 'ADDITIONAL phosphate groups on the inositol ring beyond the '
                 'one connecting to the glycerol backbone\n'
                 '- Most false positives are PI molecules with only one '
                 'phosphate group (the one connecting to glycerol)\n'
                 '\n'
                 '2. False Negatives Analysis:\n'
                 '- The code missed D-myo-Ins(1,4,5)P3 '
                 'hexakis(butyryloxymethyl) ester because it was looking too '
                 'strictly for a glycerol backbone\n'
                 '- Some phosphoinositides may have modified structures that '
                 "don't match the exact glycerol-phosphate pattern\n"
                 '\n'
                 '3. Improvements needed:\n'
                 '\n'
                 'a) Better phosphate counting:\n'
                 '- Need to distinguish between the phosphate connecting to '
                 'glycerol and additional phosphates\n'
                 '- Should count only phosphates that are ON the inositol ring '
                 'hydroxyls\n'
                 '- Should exclude the bridging phosphate from the count\n'
                 '\n'
                 'b) More flexible structure matching:\n'
                 '- Allow for modified structures while maintaining core '
                 'requirements\n'
                 '- Better handling of variations in connectivity\n'
                 '\n'
                 "Here's my proposed improvement:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 27,
    'num_false_positives': 100,
    'num_true_negatives': 104590,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.2125984251968504,
    'recall': 1.0,
    'f1': 0.35064935064935066,
    'accuracy': 0.999045045217109,
    'negative_predictive_value': 1.0}