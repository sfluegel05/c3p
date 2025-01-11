"""
Classifies: CHEBI:26935 tetraterpenoid
"""
"""
Classifies: CHEBI:26764 tetraterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.
    Tetraterpenoids are derived from tetraterpenes (C40 skeleton) and may have 
    modifications like rearrangements or removal of methyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetraterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbons - expanded range to include glycosylated derivatives
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 25 or c_count > 75:
        return False, f"Carbon count ({c_count}) outside typical range for tetraterpenoids (25-75)"

    # Count methyl groups - tetraterpenoids typically have multiple
    methyl_pattern = Chem.MolFromSmarts("[CH3]")
    methyl_count = len(mol.GetSubstructMatches(methyl_pattern))
    if methyl_count < 4:
        return False, f"Too few methyl groups ({methyl_count}) for tetraterpenoid"

    # Check for characteristic conjugated polyene patterns
    polyene_patterns = [
        "C=CC=CC=CC=C",  # Long conjugated chain
        "C=CC=CC=CC=CC=C",  # Extended conjugation
        "C(C)(C)=CC=C"  # Typical end group
    ]
    
    found_polyene = False
    for pattern in polyene_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_polyene = True
            break
    
    if not found_polyene:
        return False, "Missing characteristic polyene system"

    # Count double bonds - tetraterpenoids typically have many
    double_bond_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("C=C")))
    if double_bond_count < 8:
        return False, f"Too few double bonds ({double_bond_count}) for tetraterpenoid"

    # Check molecular weight - adjusted range for glycosylated forms
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 1200:
        return False, f"Molecular weight ({mol_wt}) outside typical range for tetraterpenoids"

    # Count rings
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count > 6:
        return False, f"Too many rings ({ring_count}) for typical tetraterpenoid"
    
    # Count nitrogens - tetraterpenoids rarely contain nitrogen
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count > 0:
        return False, f"Contains nitrogen, unusual for tetraterpenoid"
    
    # Count oxygens - increased limit for glycosylated forms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count > 15:
        return False, f"Too many oxygens ({o_count}) for typical tetraterpenoid"

    # Check for branching - more flexible pattern
    branching_patterns = [
        "[*]([*])([*])[*]",  # General branching
        "C([C,O])([C,O])[C,O]",  # Carbon branching
        "C(=C)([C,O])[C,O]"  # Branching at double bonds
    ]
    
    found_branching = False
    for pattern in branching_patterns:
        if len(mol.GetSubstructMatches(Chem.MolFromSmarts(pattern))) > 0:
            found_branching = True
            break
            
    if not found_branching:
        return False, "Insufficient characteristic branching"

    # Calculate degree of unsaturation
    du = rdMolDescriptors.CalcNumRotatableBonds(mol) + ring_count + double_bond_count
    if du < 8:
        return False, f"Insufficient degree of unsaturation ({du}) for tetraterpenoid"

    # Look for characteristic end groups
    end_groups = [
        "CC(C)=C",  # Typical isoprene end
        "C1C(C)=CCCC1(C)C",  # Cyclic end group
        "CC(=O)C",  # Keto end group
        "CC(O)=C"  # Hydroxy end group
    ]
    
    found_end_group = False
    for pattern in end_groups:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_end_group = True
            break
            
    if not found_end_group:
        return False, "Missing characteristic end groups"

    ring_info = f" and {ring_count} rings" if ring_count > 0 else ""
    oxygen_info = f" Contains {o_count} oxygen atoms." if o_count > 0 else ""

    return True, (f"Matches tetraterpenoid pattern with {c_count} carbons, "
                 f"{methyl_count} methyl groups, {double_bond_count} double bonds{ring_info}.{oxygen_info}")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26935',
                          'name': 'tetraterpenoid',
                          'definition': 'Any terpenoid derived from a '
                                        'tetraterpene. The term includes '
                                        'compounds in which the C40 skeleton '
                                        'of the parent tetraterpene has been '
                                        'rearranged or modified by the removal '
                                        'of one or more skeletal atoms '
                                        '(generally methyl groups).',
                          'parents': ['CHEBI:26873'],
                          'xrefs': ['LIPID_MAPS_class:LMPR0107'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.7203389830508474 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'C\\C(\\C=C\\C=C(/C)\\C=C\\C1=C(C)C[C@@H](O)CC1(C)C)=C/C=C/C=C(\\C)/C=C/C=C(\\C)/C=C/C1=C(C)C(=O)CCC1(C)C '
               "NAME: 3'-hydroxyechinenone REASON: CORRECT Matches "
               'tetraterpenoid pattern with 40 carbons, 10 methyl groups, 11 '
               'double bonds and 2 rings. Contains 2 oxygen atoms.\n'
               ' * SMILES: '
               'O1[C@@]2(C(C[C@]1(C[C@]2(O)C)[H])(C)C)/C=C/C(/C)=C/C=C/C(/C)=C/C=C/C=C(/C=C/C=C(/C=C/C=3C(C[C@H](O)CC3C)(C)C)\\C)\\C '
               'NAME: Cucurbitaxanthin A REASON: CORRECT Matches '
               'tetraterpenoid pattern with 40 carbons, 10 methyl groups, 10 '
               'double bonds and 3 rings. Contains 3 oxygen atoms.\n'
               ' * SMILES: '
               'CC(C)=CCC\\C(C)=C\\C=C\\C(C)=C\\C=C\\C(C)=C\\C=C\\C=C(C)\\C=C\\C=C(C)\\C=C\\C1=C(C)CCCC1(C)C '
               'NAME: gamma-carotene REASON: CORRECT Matches tetraterpenoid '
               'pattern with 40 carbons, 10 methyl groups, 12 double bonds and '
               '1 rings.\n'
               ' * SMILES: '
               'O=C(O)/C(=C/C=C/C(=C/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C1=C(CC[C@H](C1(C)C)O)C)\\C)\\C)/C)/C)/C)/C '
               'NAME: 2-hydroxytorularhodin REASON: CORRECT Matches '
               'tetraterpenoid pattern with 40 carbons, 9 methyl groups, 13 '
               'double bonds and 1 rings. Contains 3 oxygen atoms.\n'
               ' * SMILES: '
               'O=C1C(=C(/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C=C(/C=C/C(O)C(O[C@@H]2OC[C@@H](O)[C@@H]([C@H]2O)O)(C)C)\\C)\\C)\\C)/C)/C)C(C)(C)CC1O)C '
               "NAME: 1'-beta-xylosyl-3,2'-dihydroxy-beta,psi-caroten-4-one "
               'REASON: CORRECT Matches tetraterpenoid pattern with 45 '
               'carbons, 10 methyl groups, 12 double bonds and 2 rings. '
               'Contains 8 oxygen atoms.\n'
               ' * SMILES: '
               'CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\C=C\\C(C)=C\\C=C\\C=C(C)\\C=C\\C=C(C)\\C=C\\C1C(C)=CCCC1(C)C '
               'NAME: alpha-zeacarotene REASON: CORRECT Matches tetraterpenoid '
               'pattern with 40 carbons, 10 methyl groups, 11 double bonds and '
               '1 rings.\n'
               ' * SMILES: '
               'CC(\\C=C\\C=C(C)\\C=C\\C1C(C)=CC(O)CC1(C)C)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/C1C(C)=CC(O)CC1(C)C '
               'NAME: tunaxanthin REASON: CORRECT Matches tetraterpenoid '
               'pattern with 40 carbons, 10 methyl groups, 11 double bonds and '
               '2 rings. Contains 2 oxygen atoms.\n'
               ' * SMILES: '
               'O[C@@H]1CC(C(=C(C1)C)/C=C/C(/C)=C/C=C/C(/C)=C/C=C/C=C(/C=C/C=C(/C=C/C=2C(C[C@@H](O)CC2C)(C)C)\\C)\\C)(C)C '
               "NAME: (3S,3'S,all-E)-Zeaxanthin REASON: CORRECT Matches "
               'tetraterpenoid pattern with 40 carbons, 10 methyl groups, 11 '
               'double bonds and 2 rings. Contains 2 oxygen atoms.\n'
               ' * SMILES: '
               'O1C(C1CC/C(/C)=C\\C=C\\C(\\C)=C\\C=C\\C(=C/C=C/C=C(/C=C/C=C(\\C=C\\C2C(CCC=C2C)(C)C)/C)\\C)\\C)(C)C '
               'NAME: delta-Carotene-1,2-epoxide REASON: CORRECT Matches '
               'tetraterpenoid pattern with 40 carbons, 10 methyl groups, 11 '
               'double bonds and 2 rings. Contains 1 oxygen atoms.\n'
               ' * SMILES: '
               'O=C1C(=C(C=CC(=CC=CC(=CC=CC=C(C=CC=C(C=CC2=C(CC(O)CC2(C)C)C)C)C)C)C)C(C)(C)CC1C3OC(C(O)C(C3O)O)CO)C '
               "NAME: (3S, 3'R)-adonixanthin-beta-D-glucoside REASON: CORRECT "
               'Matches tetraterpenoid pattern with 46 carbons, 10 methyl '
               'groups, 11 double bonds and 3 rings. Contains 7 oxygen atoms.\n'
               ' * SMILES: '
               'O[C@@]1(C(C[C@H](O)C[C@]1(O)C)(C)C)/C=C/C(/C)=C/C=C/C(/C)=C/C=C/C=C(/C=C/C=C(/C=C/[C@@]2(O)C(C[C@H](O)C[C@]2(O)C)(C)C)\\C)\\C '
               'NAME: Mactraxanthin REASON: CORRECT Matches tetraterpenoid '
               'pattern with 40 carbons, 10 methyl groups, 9 double bonds and '
               '2 rings. Contains 6 oxygen atoms.\n'
               ' * SMILES: '
               'O=C/C(=C/C=C/C(=C/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C1=C(CCCC1(C)C)C)\\C)\\C)/C)/C)/C)/C '
               'NAME: Torularhodinaldehyde REASON: CORRECT Matches '
               'tetraterpenoid pattern with 40 carbons, 9 methyl groups, 13 '
               'double bonds and 1 rings. Contains 1 oxygen atoms.\n'
               ' * SMILES: '
               'COC(C)(C)C\\C=C\\C(C)=C\\C=C\\C(C)=C\\C=C\\C(C)=C\\C=C\\C=C(C)\\C=C\\C=C(C)\\C=C\\C=C(/C)CCC=C(C)C '
               'NAME: Anhydrorhodovibrin REASON: CORRECT Matches '
               'tetraterpenoid pattern with 41 carbons, 11 methyl groups, 13 '
               'double bonds. Contains 1 oxygen atoms.\n'
               ' * SMILES: '
               'S(O[C@H]1CC(C(=C(C1=O)C)/C=C/C(/C)=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C=2C([C@@H](O)[C@H](O)CC2C)(C)C)\\C)\\C)/C)(C)C)(O)(=O)=O '
               'NAME: Erythroxanthin sulfate REASON: CORRECT Matches '
               'tetraterpenoid pattern with 40 carbons, 10 methyl groups, 11 '
               'double bonds and 2 rings. Contains 7 oxygen atoms.\n'
               ' * SMILES: '
               'OC1CC([C@@H](C(=C1)C)/C=C/C(/C)=C/C=C/C(/C)=C/C=C/C=C(/C=C/C=C(\\C=C\\C=2C(C[C@H](O)CC2C)(C)C)/C)\\C)(C)C '
               "NAME: (3R,3'R,6'R,9-cis)-Carotene-3,3'-diol REASON: CORRECT "
               'Matches tetraterpenoid pattern with 40 carbons, 10 methyl '
               'groups, 11 double bonds and 2 rings. Contains 2 oxygen atoms.\n'
               ' * SMILES: '
               'S(O[C@H]1CC(C(=C(C1)C)/C=C/C(/C)=C/C=C/C(/C)=C/C=C/C=C(/C=C/C=C(/C=C/C=2C([C@@H](O)[C@H](O)CC2C)(C)C)\\C)\\C)(C)C)(O)(=O)=O '
               'NAME: Caloxanthin sulfate REASON: CORRECT Matches '
               'tetraterpenoid pattern with 40 carbons, 10 methyl groups, 11 '
               'double bonds and 2 rings. Contains 6 oxygen atoms.\n'
               ' * SMILES: '
               'CC(\\C=C\\C=C(C)\\C=C\\C(=O)[C@]1(C)C[C@@H](O)CC1(C)C)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/C(=O)[C@]1(C)C[C@@H](O)CC1(C)C '
               'NAME: capsorubin REASON: CORRECT Matches tetraterpenoid '
               'pattern with 40 carbons, 10 methyl groups, 9 double bonds and '
               '2 rings. Contains 4 oxygen atoms.\n'
               ' * SMILES: '
               'CO[C@@H]1[C@H](C)O[C@@H](O[C@@H](\\C=C\\C(C)=C\\C=C\\C(C)=C\\C=C\\C(C)=C\\C=C\\C=C(C)\\C=C\\C=C(C)\\C=C\\C2=C(C)CCCC2(C)C)C(C)(C)O)[C@@H](OC)[C@@H]1O '
               "NAME: (2'S)-Deoxymyxol 2'-(2,4-di-O-methyl-alpha-L-fucoside) "
               'REASON: CORRECT Matches tetraterpenoid pattern with 48 '
               'carbons, 13 methyl groups, 12 double bonds and 2 rings. '
               'Contains 6 oxygen atoms.\n'
               ' * SMILES: '
               'CC(\\C=C\\C=C(C)\\C=C\\[C@@H]1C(C)=CCCC1(C)C)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/C1=C(C)CCCC1(C)C '
               "NAME: (6'S)-beta,epsilon-carotene REASON: CORRECT Matches "
               'tetraterpenoid pattern with 40 carbons, 10 methyl groups, 11 '
               'double bonds and 2 rings.\n'
               ' * SMILES: '
               'CC(C)=CCC\\C(C)=C\\C=C\\C(C)=C\\C=C\\C(C)=C\\C=C\\C=C(CO)/C=C/C=C(C)/C=C/C=C(\\C)CCCC(C)(C)O '
               'NAME: rhodopinol REASON: CORRECT Matches tetraterpenoid '
               'pattern with 40 carbons, 9 methyl groups, 12 double bonds. '
               'Contains 2 oxygen atoms.\n'
               ' * SMILES: '
               'CC(C)=CCC\\C(C)=C\\C=C\\C(C)=C\\C=C\\C(C)=C\\C=C\\C=C(C)\\C=C\\C=C(C)\\C=C\\C1=C(C)C[C@@H](O)CC1(C)C '
               'NAME: rubixanthin REASON: CORRECT Matches tetraterpenoid '
               'pattern with 40 carbons, 10 methyl groups, 12 double bonds and '
               '1 rings. Contains 1 oxygen atoms.\n'
               ' * SMILES: '
               'CC(\\C=C\\C=C(C)\\C=C\\[C@@]12O[C@]1(C)C[C@@H](O)CC2(C)C)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/[C@@]12O[C@]1(C)C[C@@H](O)CC2(C)C '
               'NAME: all-trans-violaxanthin REASON: CORRECT Matches '
               'tetraterpenoid pattern with 40 carbons, 10 methyl groups, 9 '
               'double bonds and 4 rings. Contains 4 oxygen atoms.\n'
               ' * SMILES: '
               'OC1CC(C(=C(C1)C)/C=C/C(/C)=C\\C=C\\C(\\C)=C\\C=C/C=C(/C=C\\C=C(/C=C/C=C(/CCC=C(C)C)\\C)\\C)\\C)(C)C '
               'NAME: (all-E)-Rubixanthin REASON: CORRECT Matches '
               'tetraterpenoid pattern with 40 carbons, 10 methyl groups, 12 '
               'double bonds and 1 rings. Contains 1 oxygen atoms.\n'
               ' * SMILES: '
               'O=C/C(/C=C/C=C(/C=C/C=C(/C=C/CC(OC)(C)C)\\C)\\C)=C/C=C/C=C(/C=C/C=C(/C=C/C=C(/C=C/CC(OC)(C)C)\\C)\\C)\\C '
               'NAME: Rhodobacterioxanthin REASON: CORRECT Matches '
               'tetraterpenoid pattern with 42 carbons, 11 methyl groups, 13 '
               'double bonds. Contains 3 oxygen atoms.\n'
               ' * SMILES: '
               'OC/C(=C/CC1C(C(C(C)=CC1)/C=C/C(=C/C=C/C(=C/C=C/C=C(/CC/C=C(/CC/C=C(/CCC=C(C)C)\\C)\\C)\\C)/C)/C)(C)C)/C '
               'NAME: Nonaprenoxanthin REASON: CORRECT Matches tetraterpenoid '
               'pattern with 45 carbons, 11 methyl groups, 11 double bonds and '
               '1 rings. Contains 1 oxygen atoms.\n'
               'False positives: SMILES: '
               'CC(\\C=C\\C1=C(C)CCCC1(C)C)=C/C=C/C(C)=C/c1cc(\\C=C\\C=C(C)\\C=C\\C2=C(C)CCCC2(C)C)cc[n+]1CCO '
               'NAME: N-retinylidene-N-retinylethanolamine REASON: WRONGLY '
               'CLASSIFIED Matches tetraterpenoid pattern with 42 carbons, 9 '
               'methyl groups, 9 double bonds and 3 rings. Contains 1 oxygen '
               'atoms.\n'
               ' * SMILES: C1(=C(C(CCC1(C)C)=O)C)/C=C/C(=C/C=C/C(=C/CO)/C)/C '
               'NAME: all-trans-4-oxoretinol REASON: WRONGLY CLASSIFIED '
               'Matches tetraterpenoid pattern with 20 carbons, 5 methyl '
               'groups, 5 double bonds and 1 rings. Contains 2 oxygen atoms.\n'
               ' * SMILES: '
               'O=C1C=C(/C(=C/C=C(/C=C/C=C(/C=C/C=C/C(=C/C=C/C(=O)C)/C)\\C)\\C)/C(C1)(C)C)C '
               "NAME: (6E)-9'-apo-rhodoxanthinone REASON: WRONGLY CLASSIFIED "
               'Matches tetraterpenoid pattern with 29 carbons, 7 methyl '
               'groups, 9 double bonds and 1 rings. Contains 2 oxygen atoms.\n'
               ' * SMILES: '
               'C1(C)(C)CCC(C(=C1\\C=C\\C(=C\\C=C\\C(=C\\C(=O)[O-])\\C)\\C)CO)O '
               'NAME: all-trans-4,18-dihydroxyretinoate REASON: WRONGLY '
               'CLASSIFIED Matches tetraterpenoid pattern with 20 carbons, 4 '
               'methyl groups, 5 double bonds and 1 rings. Contains 4 oxygen '
               'atoms.\n'
               ' * SMILES: '
               'O[C@H]1CC(C(=C(C1)C)\\C=C\\C(\\C)=C\\C=C\\C(\\C)=C\\C=C\\C=C(\\C=C\\C=C(\\C#CC=2C(C[C@H](O)CC2C)(C)C)/C)/C)(C)C '
               'NAME: Diatoxanthin REASON: WRONGLY CLASSIFIED Matches '
               'tetraterpenoid pattern with 40 carbons, 10 methyl groups, 10 '
               'double bonds and 2 rings. Contains 2 oxygen atoms.\n'
               ' * SMILES: '
               'O1[C@@H]2[C@@](O)(C[C@H](NC(=O)/C=C/C=C/C=C/[C@@H](C[C@@H](C/C(=C/[C@@H](CC)C)/C)C)C)CO)C(=O)C=C[C@H]12 '
               'NAME: (+)-scyphostatin REASON: WRONGLY CLASSIFIED Matches '
               'tetraterpenoid pattern with 29 carbons, 5 methyl groups, 5 '
               'double bonds and 2 rings. Contains 5 oxygen atoms.\n'
               ' * SMILES: '
               'C1(C)(C)C(\\C=C\\C(=C\\C=C\\C(=C\\COC(CCCCCC)=O)\\C)\\C)=C(C)CCC1 '
               'NAME: all-trans-retinyl heptanoate REASON: WRONGLY CLASSIFIED '
               'Matches tetraterpenoid pattern with 27 carbons, 6 methyl '
               'groups, 5 double bonds and 1 rings. Contains 2 oxygen atoms.\n'
               ' * SMILES: '
               'O=C1O[C@H]2C(C)C(C1([C@@]2(/C=C/C=C/C=C/C=C/C=C(/[C@H]3O[C@@H](C(C)=C[C@]3(O)C)C)\\C)C)C)=O '
               'NAME: Prugosene A2 REASON: WRONGLY CLASSIFIED Matches '
               'tetraterpenoid pattern with 28 carbons, 7 methyl groups, 6 '
               'double bonds and 3 rings. Contains 5 oxygen atoms.\n'
               ' * SMILES: '
               'O=C1O[C@H]([C@H]([C@@H](O)C2=CC(O)=CC=C2)C)C[C@@H](O)CC[C@H](C=CC=C(C)[C@@H](CC=CC=CC=CC=C1C)O)C '
               'NAME: Antarlide E REASON: WRONGLY CLASSIFIED Matches '
               'tetraterpenoid pattern with 33 carbons, 4 methyl groups, 6 '
               'double bonds and 2 rings. Contains 6 oxygen atoms.\n'
               ' * SMILES: '
               'O=C(O[C@H]1[C@@H](OC)[C@@H](C(=C)CC1)C2(O[C@@H]2CC=C(C)C)C)/C=C/C=C/C=C/C=C/C(=O)O '
               'NAME: 2,4,6,8-Decatetraenedioic acid, '
               'mono[2-methoxy-4-methylene-3-[2-methyl-3-(3-methyl-2-butenyl)-2-oxiranyl]cyclohex-1-yl] '
               'ester [1R-[1beta(2E,4E,6E,8E),2beta,3alpha(2R∗,3R∗)]] (2) '
               'REASON: WRONGLY CLASSIFIED Matches tetraterpenoid pattern with '
               '26 carbons, 4 methyl groups, 6 double bonds and 2 rings. '
               'Contains 6 oxygen atoms.\n'
               ' * SMILES: '
               'C1(C)(CO)CCC(C(=C1\\C=C\\C(=C\\C=C\\C(=C\\C(=O)O)\\C)\\C)C)=O '
               'NAME: all-trans-4-oxo-16-hydroxyretinoic acid REASON: WRONGLY '
               'CLASSIFIED Matches tetraterpenoid pattern with 20 carbons, 4 '
               'methyl groups, 5 double bonds and 1 rings. Contains 4 oxygen '
               'atoms.\n'
               ' * SMILES: '
               'OC1CC(C(C(=C1)C)/C=C/C(/C)=C/C=C/C(=C/C=C/C=C(/CCCO)\\C)/C)(C)C '
               'NAME: alpha-Micropteroxanthin B REASON: WRONGLY CLASSIFIED '
               'Matches tetraterpenoid pattern with 27 carbons, 6 methyl '
               'groups, 7 double bonds and 1 rings. Contains 2 oxygen atoms.\n'
               ' * SMILES: '
               'ClC1=C2C(=COC(=C2)/C=C/[C@H](CC)C)[C@@H]3C=4C(=O)O[C@@H]([C@H](C4O[C@@]3(C1=O)C)C)C '
               'NAME: Chaetomugilin F REASON: WRONGLY CLASSIFIED Matches '
               'tetraterpenoid pattern with 23 carbons, 5 methyl groups, 5 '
               'double bonds and 4 rings. Contains 5 oxygen atoms.\n'
               ' * SMILES: '
               'O=C(O)/C=C/C=C/C=C/C=C/[C@@]1([C@@H](C(=C[C@H]2[C@H]1CC[C@@H]([C@@H]2O)C)CO)/C(=C/C)/C)C '
               'NAME: Janthinopolyenemycin B REASON: WRONGLY CLASSIFIED '
               'Matches tetraterpenoid pattern with 26 carbons, 4 methyl '
               'groups, 6 double bonds and 2 rings. Contains 4 oxygen atoms.\n'
               ' * SMILES: '
               'O1[C@@H](OC([C@@H](O)/C=C/C=C/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C2=C(C(=CC=C2C)C)C)\\C)\\C)/C)/C)(C)C)[C@H](O)[C@@H](O)[C@@H]([C@H]1CO)O '
               'NAME: Sioxanthin REASON: WRONGLY CLASSIFIED Matches '
               'tetraterpenoid pattern with 45 carbons, 9 methyl groups, 11 '
               'double bonds and 2 rings. Contains 7 oxygen atoms.\n'
               ' * SMILES: '
               'O=C1C=C[C@H]2[C@@H]([C@@]1(O)C[C@H](NC(=O)/C=C/C=C/C=C/C(CC(C/C(=C/C(CC)C)/C)C)C)CO)O2 '
               'NAME: Scyphostatin REASON: WRONGLY CLASSIFIED Matches '
               'tetraterpenoid pattern with 29 carbons, 5 methyl groups, 5 '
               'double bonds and 2 rings. Contains 5 oxygen atoms.\n'
               ' * SMILES: '
               'O=C(O)/C=C/C=C/C(=C\\C=C\\C1[C@](O)(/C(=C/C)/C)[C@H]2O[C@H]2[C@H]3[C@H]1CC=C(C3)C)/C '
               'NAME: Anthracobic acid B REASON: WRONGLY CLASSIFIED Matches '
               'tetraterpenoid pattern with 25 carbons, 4 methyl groups, 6 '
               'double bonds and 3 rings. Contains 4 oxygen atoms.\n'
               ' * SMILES: '
               'O=C1C2=C(C3=CC=4C(=C(C)C(=C(C4O[C@]3([C@@]1(O)C)O)O)O)CO)C=C(/C=C/C(=C/[C@H](CC)C)/C)OC2 '
               'NAME: Azasperpyranone A REASON: WRONGLY CLASSIFIED Matches '
               'tetraterpenoid pattern with 28 carbons, 5 methyl groups, 5 '
               'double bonds and 4 rings. Contains 8 oxygen atoms.\n'
               ' * SMILES: '
               'OC1CC(C(=C(C1)C)/C=C/C(/C)=C\\C=C\\C(\\C)=C\\C=C/C=C(/C=C/C=C(\\C)/C(=O)C)\\C)(C)C '
               'NAME: 3-Hydroxysintaxanthin REASON: WRONGLY CLASSIFIED Matches '
               'tetraterpenoid pattern with 31 carbons, 8 methyl groups, 9 '
               'double bonds and 1 rings. Contains 2 oxygen atoms.\n'
               ' * SMILES: '
               'O1[C@@H](OC(CCC/C(/C)=C/C=C/C(/C)=C/C=C/C=C(\\C)/C=C/C=C(/C=C/C=C(C)C)\\C)(C)C)C(O)C(O)[C@H](O)[C@H]1CO '
               'NAME: OH-Diaponeurosporene glucoside ester REASON: WRONGLY '
               'CLASSIFIED Matches tetraterpenoid pattern with 36 carbons, 8 '
               'methyl groups, 9 double bonds and 1 rings. Contains 6 oxygen '
               'atoms.\n'
               ' * SMILES: '
               'O=C1C=C(C)[C@@H]2[C@H]1C(=CC[C@H]3/C(=C(/C=C/C=C(C)C)\\C)/CC[C@@]3(C2)C)C=O '
               'NAME: 14,15-dehydro-6-epi-ophiobolin G REASON: WRONGLY '
               'CLASSIFIED Matches tetraterpenoid pattern with 25 carbons, 5 '
               'methyl groups, 5 double bonds and 3 rings. Contains 2 oxygen '
               'atoms.\n'
               ' * SMILES: '
               'COC1CC(O)\\C=C\\C=C\\C=C/C(C)C(OC(=O)\\C=C/C=C/C=C\\C(C)=C\\C(C)C(O)\\C=C/C=C/c2coc(n2)C1C)C(C)C(O)CC(C)O '
               'NAME: Chivosazole F REASON: WRONGLY CLASSIFIED Matches '
               'tetraterpenoid pattern with 41 carbons, 7 methyl groups, 9 '
               'double bonds and 2 rings. Contains 8 oxygen atoms.\n'
               ' * SMILES: '
               'O=C(O)/C=C/C=C/C(=C/C=C/C1[C@](O)(/C(=C/C)/C)[C@H]2O[C@H]2[C@H]3[C@H]1CC=C(C3)C)/C '
               'NAME: Anthracobic acid A REASON: WRONGLY CLASSIFIED Matches '
               'tetraterpenoid pattern with 25 carbons, 4 methyl groups, 6 '
               'double bonds and 3 rings. Contains 4 oxygen atoms.\n'
               ' * SMILES: '
               'C\\C(CO)=C/C=C/C(C)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/C1=C(C)C[C@@H](O)CC1(C)C '
               "NAME: (3R)-3-hydroxy-8'-apo-beta-carotenol REASON: WRONGLY "
               'CLASSIFIED Matches tetraterpenoid pattern with 30 carbons, 7 '
               'methyl groups, 9 double bonds and 1 rings. Contains 2 oxygen '
               'atoms.\n'
               ' * SMILES: CC(\\C=C\\C=C(C)\\C=C\\C1=C(C)CC(O)CC1(C)C)=C/C=O '
               'NAME: all-trans-3-Hydroxyretinal REASON: WRONGLY CLASSIFIED '
               'Matches tetraterpenoid pattern with 20 carbons, 5 methyl '
               'groups, 5 double bonds and 1 rings. Contains 2 oxygen atoms.\n'
               'False negatives: SMILES: '
               'CC(\\C=C\\C=C(/C)C(O)=O)=C/C=C/C=C(C)/C=C/C=C(\\C)C(O)=O NAME: '
               'crocetin REASON: MISSED Insufficient branching for '
               'tetraterpenoid\n'
               ' * SMILES: '
               'C(/C=C/C=C(/C=C/C=C(C)C)\\C)(=C\\C=C\\C(=C\\C=C\\C=C(\\C=C\\C=C(\\C=C\\C=C(\\C=C\\C=C(C)C)/C)/C)/C)\\C)\\C '
               'NAME: Bisdehydrolycopene/ Tetradehydrolycopene REASON: MISSED '
               'Insufficient branching for tetraterpenoid\n'
               ' * SMILES: '
               'O=C1C(=C(/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C2=C(C(=O)[C@@H](O)CC2(C)C)C)\\C)\\C)/C)/C)C(C)(C)C[C@@H]1OC3OC(C(O)C(C3O)O)CO)C '
               "NAME: (3S,3'S)-Astaxanthin-beta-D-glucoside REASON: MISSED Too "
               'many oxygens (9) for typical tetraterpenoid\n'
               ' * SMILES: '
               'CC(\\C=C\\C=C(C)\\C=C\\C1=C(C)C[C@H](CC1(C)C)O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/C1=C(C)C[C@H](CC1(C)C)O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O '
               'NAME: zeaxanthin bis(beta-D-glucoside) REASON: MISSED Carbon '
               'count (52) outside typical range for tetraterpenoids (20-50)\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCC(=O)O[C@@H]1CC(C)=C(\\C=C\\C(C)=C\\C=C\\C(C)=C\\C=C\\C=C(C)\\C=C\\C=C(C)\\C=C\\[C@H]2C(C)=C[C@@H](CC2(C)C)OC(=O)CCCCCCCCCCCCCCC)C(C)(C)C1 '
               'NAME: Helenien REASON: MISSED Carbon count (72) outside '
               'typical range for tetraterpenoids (20-50)\n'
               ' * SMILES: '
               'O([C@H]1CC(C(=C(C1)C)/C=C/C(/C)=C/C=C/C(/C)=C/C=C/C=C(/C=C/C=C(/C=C/C=2C(C[C@H](O[C@@H]3O[C@@H]([C@@H](O)C(O)C3O)CO)CC2C)(C)C)\\C)\\C)(C)C)[C@@H]4O[C@@H]([C@@H](O)C(O)C4O)CO '
               'NAME: Zeaxanthin diglucoside/ Zeaxanthin beta-D-diglucoside '
               'REASON: MISSED Carbon count (52) outside typical range for '
               'tetraterpenoids (20-50)\n'
               ' * SMILES: '
               'O([C@H]1CC(C(=C(C1=O)C)/C=C/C(/C)=C/C=C/C(/C)=C/C=C/C=C(/C=C/C=C(/C=C/C=2C(C[C@H](O)C(=O)C2C)(C)C)\\C)\\C)(C)C)[C@@H]3OC([C@@H](O)C(O)C3O)CO '
               'NAME: Astaxanthin glucoside REASON: MISSED Too many oxygens '
               '(9) for typical tetraterpenoid\n'
               ' * SMILES: '
               'O=C(OCC1OC(O[C@@H]2CC(=C(/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C3=C(CCCC3(C)C)C)\\C)\\C)/C)/C)C(C2)(C)C)C)C(O)C(C1O)O)CCCCCCCCCC(C)C '
               'NAME: C1-13 thermocryptoxanthin-13 REASON: MISSED Carbon count '
               '(59) outside typical range for tetraterpenoids (20-50)\n'
               ' * SMILES: '
               'CCC(C)CCCCCCCCCCC(=O)OC[C@H]1O[C@@H](OC(=O)C(\\C)=C\\C=C\\C(C)=C\\C=C\\C(C)=C\\C=C\\C=C(C)\\C=C\\C=C(/C)CCC=C(C)C)[C@H](O)[C@@H](O)[C@@H]1O '
               'NAME: staphyloxanthin REASON: MISSED Carbon count (51) outside '
               'typical range for tetraterpenoids (20-50)\n'
               ' * SMILES: '
               'O=C1C(=C(/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C2=C(C(=O)[C@@H](OC3OC(C(O)C(C3O)O)C)CC2(C)C)C)\\C)\\C)/C)/C)C(C)(C)C[C@@H]1OC4OC(C(O)C(C4O)O)C)C '
               'NAME: Astaxanthin dirhamnoside REASON: MISSED Carbon count '
               '(52) outside typical range for tetraterpenoids (20-50)\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'Cc1cc(C)c2Cc3c(C)cc(C)c(Cc4c(C)cc(C)c(Cc5c(C)cc(C)c(Cc1c2C)c5C)c4C)c3C',
                                     'name': '4,6,10,12,16,18,22,24,25,26,27,28-dodecamethylcalix[4]arene',
                                     'reason': 'Missing characteristic polyene '
                                               'system'},
                                 {   'smiles': 'O=C1C(=C2C=C3[C@]([C@@H](C(C)C)[C@@H]([C@H]3O)OC(=O)C)(C)CC[C@]2(C)CC1)COC(=O)C',
                                     'name': 'Dahliane E',
                                     'reason': 'Carbon count (24) outside '
                                               'typical range for '
                                               'tetraterpenoids (25-75)'},
                                 {   'smiles': 'O[C@@H]([C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)N[C@@H](CC(C)C)C(O)=O)C',
                                     'name': 'Glu-Thr-Leu',
                                     'reason': 'Carbon count (15) outside '
                                               'typical range for '
                                               'tetraterpenoids (25-75)'},
                                 {   'smiles': 'CCOc1ccc(NC(=O)C(C)O)cc1',
                                     'name': 'p-Lactophenetide',
                                     'reason': 'Carbon count (11) outside '
                                               'typical range for '
                                               'tetraterpenoids (25-75)'},
                                 {   'smiles': 'O=C1NCC=2C1=C3C(N([C@@H]4O[C@H]([C@@H](O)[C@H]([C@H]4OC)O)C)C5=C3C=CC=C5)=C6NC7=C(C26)C=CC=C7',
                                     'name': "3'-epi-5'-methoxy-K252d",
                                     'reason': 'Too few methyl groups (2) for '
                                               'tetraterpenoid'},
                                 {   'smiles': 'O=C(OC1C(O)C(OC(C1O)CO)OCC2OC(OCC(OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)C(O)C(C2O)O)CCCCCCCCCCCCCCC',
                                     'name': '[(2S)-2-hexadecanoyloxy-3-[(2S,3R,4S,5S,6R)-6-[[(2S,3R,4S,5S,6R)-4-hexadecanoyloxy-3,5-dihydroxy-6-(hydroxymethyl)tetrahydropyran-2-yl]oxymethyl]-3,4,5-trihydroxy-tetrahydropyran-2-yl]oxy-propyl] '
                                             '(9E,12E,15E)-octadeca-9,12,15-trienoate',
                                     'reason': 'Too few methyl groups (3) for '
                                               'tetraterpenoid'},
                                 {   'smiles': 'O=C1C2=C(O)C=C(OC)C=C2C(=O)C3=C1[C@@H]([C@@H](O)[C@]([C@@H]3O)(O)C)C',
                                     'name': 'Altersolanol G',
                                     'reason': 'Carbon count (17) outside '
                                               'typical range for '
                                               'tetraterpenoids (25-75)'},
                                 {   'smiles': '[H][C@]1(O[C@](O)(C[C@H](O)[C@H]1NC(=O)CO)C([O-])=O)[C@H](O)[C@H](O)CO',
                                     'name': 'N-glycoloyl-alpha-neuraminate',
                                     'reason': 'Carbon count (11) outside '
                                               'typical range for '
                                               'tetraterpenoids (25-75)'},
                                 {   'smiles': 'OC(C(O)C/C=C\\C/C=C\\C/C=C\\CC)C/C=C\\C/C=C\\CCC(O)=O',
                                     'name': '10,11-DiHDPE',
                                     'reason': 'Carbon count (22) outside '
                                               'typical range for '
                                               'tetraterpenoids (25-75)'},
                                 {   'smiles': '[Na+].[H][C@]12SCC(C)=C(N1C(=O)[C@H]2NC(=O)[C@H](N)c1ccccc1)C([O-])=O',
                                     'name': 'cephalexin sodium',
                                     'reason': 'Carbon count (16) outside '
                                               'typical range for '
                                               'tetraterpenoids (25-75)'}],
    'sample_false_negatives': [   {   'smiles': 'CC(\\C=C\\C=C(/C)C(O)=O)=C/C=C/C=C(C)/C=C/C=C(\\C)C(O)=O',
                                      'name': 'crocetin',
                                      'reason': 'Carbon count (20) outside '
                                                'typical range for '
                                                'tetraterpenoids (25-75)'},
                                  {   'smiles': 'C\\C(C=O)=C/C=C/C=C(\\C)/C=C/C=C(\\C)C=[C@]=C1C(C)(C)C[C@H](O)C[C@@]1(C)O',
                                      'name': "(3S,5R,6R)-3,5-dihydroxy-6,7-didehydro-5,6-dihydro-12'-apo-beta-caroten-12'-al",
                                      'reason': 'Too few double bonds (7) for '
                                                'tetraterpenoid'},
                                  {   'smiles': 'C\\C(C=O)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/[C@@]12O[C@]1(C)C[C@@H](O)CC2(C)C',
                                      'name': "(3S,5R,6S)-5,6-epoxy-3-hydroxy-5,6-dihydro-12'-apo-beta-caroten-12'-al",
                                      'reason': 'Too few double bonds (6) for '
                                                'tetraterpenoid'},
                                  {   'smiles': 'O=C1C(=C(/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C=C(/CCCC(OC2OC(C(O)C(C2O)O)COC(=O)C(C(O)CCCCCCCCCCCCCCCCCCCC(/C=C\\CCCCCCCCCCCCCC(=O)O)C)CCCCCCCCCCCCCCCCCCCC)(C)C)\\C)\\C)\\C)/C)/C)C(C)(C)CC1)C',
                                      'name': "1'-[(6-O-Mycoloyl-beta-D-glucopyranosyl)",
                                      'reason': 'Carbon count (106) outside '
                                                'typical range for '
                                                'tetraterpenoids (25-75)'},
                                  {   'smiles': 'CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\C1C(COP(O)(=O)OP(O)(O)=O)C1(C)CC\\C=C(/C)CC\\C=C(/C)CCC=C(C)C',
                                      'name': 'prephytoene diphosphate',
                                      'reason': 'Missing characteristic '
                                                'polyene system'},
                                  {   'smiles': '[H]C(=O)\\C=C\\C=C(C)\\C=C\\C=C(C)\\C=C\\C1=C(C)CCCC1(C)C',
                                      'name': "14'-apo-beta-carotenal",
                                      'reason': 'Carbon count (22) outside '
                                                'typical range for '
                                                'tetraterpenoids (25-75)'},
                                  {   'smiles': 'O=C(OCC1OC(O[C@@H]2CC(=C(/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C3=C(C[C@@H](OC4C(O)C(O)C(O)C(C4)COC(=O)CCCCCCCCCC(C)C)CC3(C)C)C)\\C)\\C)/C)/C)C(C2)(C)C)C)C(O)C(C1O)O)CCCCCCCCCC(C)C',
                                      'name': 'Thermobiszeaxanthin Z2-13-13',
                                      'reason': 'Carbon count (79) outside '
                                                'typical range for '
                                                'tetraterpenoids (25-75)'},
                                  {   'smiles': 'CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\[C@@H]1[C@@H](COP(O)(=O)OP(O)(O)=O)[C@]1(C)CC\\C=C(/C)CC\\C=C(/C)CCC=C(C)C',
                                      'name': '(1R,2R,3R)-prephytoene '
                                              'diphosphate',
                                      'reason': 'Missing characteristic '
                                                'polyene system'},
                                  {   'smiles': 'CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\[C@H]1[C@H](COP(O)(=O)OP(O)(O)=O)[C@@]1(C)CC\\C=C(/C)CC\\C=C(/C)CCC=C(C)C',
                                      'name': '(1S,2S,3S)-prephytoene '
                                              'diphosphate',
                                      'reason': 'Missing characteristic '
                                                'polyene system'}],
    'attempt': 2,
    'reasoning': 'Let me analyze the issues with the previous code based on '
                 'the outcomes:\n'
                 '\n'
                 '1. False Negatives Analysis:\n'
                 '- The code missed crocetin and Bisdehydrolycopene due to '
                 '"Insufficient branching" - this suggests our branching check '
                 'is too strict\n'
                 '- Several compounds were rejected due to carbon count limits '
                 '(52-72 carbons) - we need to expand the range to include '
                 'glycosylated derivatives\n'
                 '- Some compounds were rejected due to oxygen count limit '
                 '(>8) - we need to increase this limit to accommodate '
                 'glycosylated forms\n'
                 '\n'
                 '2. False Positives Analysis:\n'
                 '- Many smaller molecules (20-30 carbons) were incorrectly '
                 'classified\n'
                 '- Several compounds with similar patterns but different core '
                 'structures were misidentified\n'
                 '- Need better structural checks for the tetraterpene core\n'
                 '\n'
                 "Here's how we can improve the code:\n"
                 '\n'
                 '1. Structural Improvements:\n'
                 '- Add pattern matching for characteristic tetraterpene '
                 'backbone fragments\n'
                 '- Better detection of carotenoid-like conjugated systems\n'
                 '- More flexible branching patterns detection\n'
                 '\n'
                 '2. Parameter Adjustments:\n'
                 '- Expand carbon count range to include glycosylated forms '
                 '(up to 75 carbons)\n'
                 '- Increase oxygen limit to accommodate glycosylated forms '
                 '(up to 15 oxygens)\n'
                 '- Adjust molecular weight range accordingly\n'
                 '\n'
                 '3. Additional Checks:\n'
                 '- Add check for characteristic end groups\n'
                 '- Improve conjugated system detection\n'
                 '- Add check for typical substitution patterns\n'
                 '\n'
                 "Here's the improved code:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 193,
    'num_false_positives': 100,
    'num_true_negatives': 108954,
    'num_false_negatives': 9,
    'num_negatives': None,
    'precision': 0.658703071672355,
    'recall': 0.9554455445544554,
    'f1': 0.7797979797979798,
    'accuracy': 0.999002343120744,
    'negative_predictive_value': 0.9999174031551995}