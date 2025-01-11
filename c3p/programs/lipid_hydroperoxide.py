"""
Classifies: CHEBI:61051 lipid hydroperoxide
"""
"""
Classifies: lipid hydroperoxide
Definition: Any lipid carrying one or more hydroperoxy substituents.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_lipid_hydroperoxide(smiles: str):
    """
    Determines if a molecule is a lipid hydroperoxide based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (is_match, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for hydroperoxy group (-OO)
    hydroperoxy_pattern = Chem.MolFromSmarts("[OX2][OX2H]")
    hydroperoxy_matches = mol.GetSubstructMatches(hydroperoxy_pattern)
    if not hydroperoxy_matches:
        return False, "No hydroperoxy (-OO) group found"

    # Check for lipid characteristics
    
    # Look for carboxylic acid or ester group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OH,OR]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid or ester group found"
    
    # Count carbons (lipids typically have >8 carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 8:
        return False, "Carbon chain too short for a lipid"
    
    # Look for carbon chain
    alkyl_chain = Chem.MolFromSmarts("[CH2][CH2][CH2]")
    if not mol.HasSubstructMatch(alkyl_chain):
        return False, "No alkyl chain found"

    # Count hydroperoxy groups
    num_hydroperoxy = len(hydroperoxy_matches)
    
    # Check for double bonds (most lipid hydroperoxides are unsaturated)
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    has_double_bonds = len(double_bond_matches) > 0
    
    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150:  # Arbitrary minimum weight for a lipid hydroperoxide
        return False, "Molecular weight too low for a lipid hydroperoxide"

    # Construct reason string
    reason_parts = []
    reason_parts.append(f"Contains {num_hydroperoxy} hydroperoxy group(s)")
    if has_double_bonds:
        reason_parts.append(f"has {len(double_bond_matches)} double bond(s)")
    reason_parts.append(f"contains {carbon_count} carbons")
    reason = ", ".join(reason_parts)
    
    return True, reason


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:61051',
                          'name': 'lipid hydroperoxide',
                          'definition': 'Any lipid carrying one or more '
                                        'hydroperoxy substituents.',
                          'parents': ['CHEBI:18059', 'CHEBI:35924'],
                          'xrefs': ['KEGG:C01025'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'C1C=C(N=C2N1NN=N2)C3=CC=CC=C3',
                                     'name': '5-phenyl-1,7-dihydrotetrazolo[1,5-a]pyrimidine',
                                     'reason': 'No hydroperoxy (-OO) group '
                                               'found'},
                                 {   'smiles': 'O(C(=O)CCCCCCCCCCC/C=C\\CCCCCCCC)[C@H](COC(=O)CCCCCCCCC/C=C\\CCCCCC)COC(=O)CCCCCCC/C=C\\CCCC',
                                     'name': 'TG(14:1(9Z)/22:1(13Z)/18:1(11Z))',
                                     'reason': 'No hydroperoxy (-OO) group '
                                               'found'},
                                 {   'smiles': 'CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@@H](O[C@@H]3O[C@H](CO[C@H]4O[C@H](CO[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H](O)[C@H](O[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H]4O)[C@@H](O)[C@H](O[C@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]4O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4NC(C)=O)[C@@H]3O)[C@H](O)[C@H]2NC(C)=O)[C@@H]1O',
                                     'name': 'beta-D-GlcpNAc-(1->2)-alpha-D-Manp-(1->3)-{alpha-D-Manp-(1->3)-[alpha-D-Manp-(1->6)]-alpha-D-Manp-(1->6)}-beta-D-Manp-(1->4)-beta-GlcpNAc-(1->4)-beta-D-GlcpNAc',
                                     'reason': 'No hydroperoxy (-OO) group '
                                               'found'},
                                 {   'smiles': '[C@@]12([C@@]([C@]([C@@H](CC1)C)(CCC(CC)C)C)(CCC[C@@H]2C)[H])C',
                                     'name': 'clerodane',
                                     'reason': 'No hydroperoxy (-OO) group '
                                               'found'},
                                 {   'smiles': 'S(O[C@@H]1[C@H](O)[C@@H](O)[C@H](O[C@@H]([C@@H](O)[C@H](O)CO[C@]2(O[C@H]([C@H](NC(=O)C)[C@@H](O)C2)[C@H](O)[C@H](O)CO)C(O)=O)[C@@H](NC(=O)C)CO)O[C@@H]1CO)(O)(=O)=O',
                                     'name': '(2R,4S,5R,6R)-5-Acetamido-2-[(2R,3S,4R,5S)-5-acetamido-4-[(2R,3R,4R,5R,6R)-3,4-dihydroxy-6-(hydroxymethyl)-5-sulfooxyoxan-2-yl]oxy-2,3,6-trihydroxyhexoxy]-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'No hydroperoxy (-OO) group '
                                               'found'},
                                 {   'smiles': 'O=C(NC(CC(C)C)C(=O)N[C@@H](CCCN=C(N)N)C=O)C(NC(=O)CC)CC(C)C',
                                     'name': 'Leupeptin Pr-LL',
                                     'reason': 'No hydroperoxy (-OO) group '
                                               'found'},
                                 {   'smiles': 'COC1=CC=C(C=C1)CCNC(=O)C(=C2C3=CC=CC=C3C(=N2)NC(=O)C4=CC=CS4)C#N',
                                     'name': 'N-[3-[1-cyano-2-[2-(4-methoxyphenyl)ethylamino]-2-oxoethylidene]-1-isoindolyl]-2-thiophenecarboxamide',
                                     'reason': 'No hydroperoxy (-OO) group '
                                               'found'},
                                 {   'smiles': 'O=C1O[C@@H](CC=2C1=C(OC)C(O)=CC2)CCC[C@@H](O)C',
                                     'name': 'Penicimarin C',
                                     'reason': 'No hydroperoxy (-OO) group '
                                               'found'},
                                 {   'smiles': 'CO[C@@H]1[C@@H]2[C@@H](C[C@@H]3[C@@](O2)([C@H]([C@@H]([C@H](O3)C(=O)OC)O)O)O)OC1N4C=NC5=C4N=CN=C5N',
                                     'name': 'LSM-4497',
                                     'reason': 'No hydroperoxy (-OO) group '
                                               'found'},
                                 {   'smiles': 'C[C@@H]1[C@H](O)CCC2=CC[C@H](C[C@]12C)C(C)=C',
                                     'name': '1-deoxycapsidiol',
                                     'reason': 'No hydroperoxy (-OO) group '
                                               'found'}],
    'sample_false_negatives': [   {   'smiles': 'C(C(O)=O)C/C=C\\C/C=C\\C/C=C\\C\\C=C/C=C/[C@@H](C/C=C\\CC)OO',
                                      'name': '(4Z,7Z,10Z,13Z,15E,17R,19Z)-17-hydroperoxydocosahexaenoic '
                                              'acid',
                                      'reason': 'No alkyl chain found'},
                                  {   'smiles': 'OC(CC/C=C\\C/C=C\\CC(\\C=C/C=C\\C=C\\[C@H](C/C=C\\CC)OO)OO)=O',
                                      'name': '(4Z,7Z,11Z,13Z,15E,17S,19Z)-10,17-bis(hydroperoxy)docosahexaenoic '
                                              'acid',
                                      'reason': 'No alkyl chain found'},
                                  {   'smiles': 'O(O)C(CCCCCCCC(OC)=O)/C=C/C(OO)/C=C/C=C\\CC',
                                      'name': 'methyl '
                                              '9,12-dihydroperoxy-10E,13E,15Z-octadecatrienoate',
                                      'reason': 'No carboxylic acid or ester '
                                                'group found'},
                                  {   'smiles': '[H]C(=CCCCCCCCCC(=O)OC)C1C2CC(OO2)C1C([H])=CC(CC)OO',
                                      'name': 'methyl '
                                              '8-[3,5-epidioxy-2-(3-hydroperoxy-1-pentenyl)-cyclopentyl]-octanoate',
                                      'reason': 'No carboxylic acid or ester '
                                                'group found'},
                                  {   'smiles': 'CC\\C=C/C\\C=C/CC(OO)\\C=C\\C=C/C\\C=C/C\\C=C/CCC(O)=O',
                                      'name': '14-HPDHE',
                                      'reason': 'No alkyl chain found'},
                                  {   'smiles': 'C(\\CC(/C=C/C=C\\C/C=C\\CCCCCO)OO)=C\\CCCC(OC)=O',
                                      'name': 'methyl '
                                              '(5Z,9E,11Z,14Z)-8-hydroperoxy-20-hydroxyicosa-5,9,11,14-tetraenoate',
                                      'reason': 'No carboxylic acid or ester '
                                                'group found'},
                                  {   'smiles': 'C(CCCCCCC[C@@H](OO)/C=C/C=C\\CCCCC)(=O)OC',
                                      'name': 'methyl '
                                              '(9R,10E,12Z)-9-hydroperoxyoctadeca-10,12-dienoate',
                                      'reason': 'No carboxylic acid or ester '
                                                'group found'},
                                  {   'smiles': 'CCCCC\\C=C/C\\C=C/C=C/[C@@H](C\\C=C/CCCC(=O)OC)OO',
                                      'name': '8(R)-HPETE methyl ester',
                                      'reason': 'No carboxylic acid or ester '
                                                'group found'},
                                  {   'smiles': 'O=C(CCC/C=C\\C/C=C\\C\\C=C/C=C/[C@H](CCCCC)OO)OC(CO)CO',
                                      'name': '15(S)-HPETE 2-glyceryl ester',
                                      'reason': 'No carboxylic acid or ester '
                                                'group found'},
                                  {   'smiles': 'C(CC/C=C\\C[C@H](OO)/C=C/C=C\\C=C\\[C@H](C/C=C\\C/C=C\\CC)OO)(O)=O',
                                      'name': '(7S,14S)-bis(hydroperoxy)-(4Z,8E,10Z,12E,16Z,19Z)-docosahexaenoic '
                                              'acid',
                                      'reason': 'No alkyl chain found'}],
    'attempt': 0,
    'reasoning': 'Let me help create a program to classify lipid '
                 'hydroperoxides. From the examples and definition, these are '
                 'lipids containing one or more -OO (hydroperoxy) groups. Key '
                 'characteristics include:\n'
                 '\n'
                 '1. Must contain at least one -OO (hydroperoxy) group\n'
                 '2. Should be a lipid (fatty acids, their derivatives, or '
                 'related compounds)\n'
                 '3. Often contains carbon chains with unsaturation (double '
                 'bonds)\n'
                 '4. Usually contains a carboxylic acid or ester group\n'
                 '\n'
                 "Here's the program:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 88,
    'num_false_positives': 11,
    'num_true_negatives': 142177,
    'num_false_negatives': 24,
    'num_negatives': None,
    'precision': 0.8888888888888888,
    'recall': 0.7857142857142857,
    'f1': 0.8341232227488152,
    'accuracy': 0.9997540407589599,
    'negative_predictive_value': 0.999831224815578}