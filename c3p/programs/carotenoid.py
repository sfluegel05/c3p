"""
Classifies: CHEBI:23044 carotenoid
"""
"""
Classifies: CHEBI:23042 carotenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_carotenoid(smiles: str):
    """
    Determines if a molecule is a carotenoid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a carotenoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get basic molecular properties
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Check for glycosylated/esterified carotenoids
    is_glycosylated = mol.HasSubstructMatch(Chem.MolFromSmarts("OC1OC(CO)C(O)C(O)C1"))
    is_esterified = mol.HasSubstructMatch(Chem.MolFromSmarts("CC(=O)OC"))
    
    # Adjust carbon count limits based on modifications
    max_carbons = 50
    if is_glycosylated:
        max_carbons = 80  # Allow for glycosylated forms
    if is_esterified:
        max_carbons = 60  # Allow for esterified forms
        
    if c_count < 20:
        return False, f"Too few carbons ({c_count}) for a carotenoid"
    if c_count > max_carbons:
        return False, f"Too many carbons ({c_count}) for a carotenoid"
    
    # Check for carotenoid backbone patterns
    backbone_patterns = [
        "C(/C=C/C(/C)=C/C=C/C(/C)=C/C=C/C=C(/C))",  # All-trans polyene
        "C(/C=C/C=C(/C)C=C/C=C(/C))",                # Shorter polyene
        "CC(C)=CCC/C(C)=C/CC/C(C)=C/CC"              # Carotenoid precursor
    ]
    
    has_backbone = False
    for pattern in backbone_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_backbone = True
            break
            
    if not has_backbone:
        return False, "No characteristic carotenoid backbone found"
    
    # Count conjugated double bonds
    conjugated_pattern = Chem.MolFromSmarts("C=CC=C")
    num_conjugated = len(mol.GetSubstructMatches(conjugated_pattern))
    
    # Look for common end groups with more specific patterns
    end_groups = [
        ("beta", "C1=C(C)C(C)(C)CCC1"),           # Beta-type end
        ("epsilon", "C1C(C)=C(C)CC1(C)C"),        # Epsilon-type end
        ("gamma", "C1=C(C)CCCC1(C)C"),            # Gamma-type end
        ("keto", "C(=O)C=CC=C"),                  # Keto end group
        ("acyclic", "C(C)(C)C=CC=C"),             # Acyclic end
        ("cyclopentyl", "C1CCCC1(C)C"),           # Cyclopentyl end
    ]
    
    found_ends = []
    for name, smarts in end_groups:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):
            found_ends.append(name)
    
    # Check for modifications
    modifications = []
    if is_glycosylated:
        modifications.append("glycosylated")
    if is_esterified:
        modifications.append("esterified")
    if mol.HasSubstructMatch(Chem.MolFromSmarts("CO")):
        modifications.append("hydroxylated")
    if mol.HasSubstructMatch(Chem.MolFromSmarts("C(=O)C")):
        modifications.append("keto")
    if mol.HasSubstructMatch(Chem.MolFromSmarts("C1OC1")):
        modifications.append("epoxidated")
    if mol.HasSubstructMatch(Chem.MolFromSmarts("C(=O)O")):
        modifications.append("carboxylated")
    
    # Count methyl branches (characteristic of carotenoids)
    methyl_pattern = Chem.MolFromSmarts("C-C(C)-C")
    num_methyls = len(mol.GetSubstructMatches(methyl_pattern))
    
    # Final classification logic
    is_likely_carotenoid = False
    reason = ""
    
    if found_ends and num_conjugated >= 4:
        is_likely_carotenoid = True
        mod_str = f" ({', '.join(modifications)})" if modifications else ""
        end_str = ", ".join(found_ends)
        reason = f"Carotenoid with {end_str} end group(s){mod_str}"
    elif num_conjugated >= 7 and num_methyls >= 4:
        is_likely_carotenoid = True
        reason = "Carotenoid based on conjugation pattern and methyl branching"
    elif c_count >= 30 and num_conjugated >= 6 and has_backbone:
        is_likely_carotenoid = True
        reason = "Carotenoid based on backbone structure"
    else:
        reason = "Does not match carotenoid structural requirements"
        
    return is_likely_carotenoid, reason


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:23044',
                          'name': 'carotenoid',
                          'definition': 'One of a class of tetraterpenoids '
                                        '(C40), formally derived from the '
                                        'acyclic parent, psi,psi-carotene by '
                                        'hydrogenation, dehydrogenation, '
                                        'cyclization, oxidation, or '
                                        'combination of these processes. This '
                                        'class includes carotenes, '
                                        'xanthophylls and certain compounds '
                                        'that arise from rearrangement of the '
                                        'skeleton of psi,psi-carotene or by '
                                        'loss of part of this structure. '
                                        'Retinoids are excluded.',
                          'parents': ['CHEBI:26935'],
                          'xrefs': ['CAS:36-88-4'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.729957805907173 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'CC(\\C=C\\C=C(C)\\C=C\\[C@H]1C(C)=C[C@H](O)CC1(C)C)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/[C@@]12O[C@]1(C)C[C@@H](O)CC2(C)C '
               'NAME: lutein 5,6-epoxide REASON: CORRECT Carotenoid with gamma '
               'end group(s) (hydroxylated, epoxidated)\n'
               ' * SMILES: '
               'CC(\\C=C\\C=C(C)\\C=C\\C1=C(C)C[C@@H](O)CC1(C)C)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/C1=C(C)C[C@@H](O)[C@H](O)C1(C)C '
               'NAME: Caloxanthin REASON: CORRECT Likely carotenoid based on '
               'carbon count and conjugation\n'
               ' * SMILES: '
               'C(=C(\\C=C\\CC(OC)(C)C)/C)/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C=C(/C=C/C(C(OC)(C)C)(O)O)\\C)\\C)\\C)/C)/C '
               'NAME: 2,2-dihydroxyspirilloxanthin REASON: CORRECT Likely '
               'carotenoid based on carbon count and conjugation\n'
               ' * SMILES: '
               'O1[C@@]2(C(C[C@]1(C[C@]2(O)C)[H])(C)C)/C=C/C(/C)=C/C=C/C(/C)=C/C=C/C=C(/C=C/C=C(/C=C/C=3C(C[C@H](O)CC3C)(C)C)\\C)\\C '
               'NAME: Cucurbitaxanthin A REASON: CORRECT Likely carotenoid '
               'based on carbon count and conjugation\n'
               ' * SMILES: '
               'OC1CC(C(=C(C1=O)C)/C=C/C(/C)=C/C=C/C(/C)=C/C=C/C=C(/C=C/C=C(/C=C/C=C(/C=C/C=C(C)C)\\C)\\C)\\C)(C)C '
               "NAME: 3-Hydroxy-3',4'-didehydro-beta,psi-caroten-4-one REASON: "
               'CORRECT Carotenoid with keto end group(s) (hydroxylated, '
               'keto)\n'
               ' * SMILES: '
               'OC1(C(C[C@@H](O)C=C1C)(C)C)/C=C/C(/C)=C/C=C/C(/C)=C/C=C/C=C(/C=C/C=C(/C(=O)CC=2C(C[C@H](O)CC2C)(C)C)\\CO)\\C '
               "NAME: 6'-Hydroxysiphonaxanthin REASON: CORRECT Carotenoid with "
               'gamma, keto end group(s) (hydroxylated, keto)\n'
               ' * SMILES: '
               'C(/C=C/C(=C/C=C/C(=C/C=C/C=C(\\C)/C=C/C=C(/C=C/[C@@]12[C@](C[C@@H](O)CC1(C)C)(C)O2)\\C)/C)/C)([C@@]3(C[C@H](CC3(C)C)O)C)=O '
               'NAME: capsanthin 5,6-epoxide REASON: CORRECT Carotenoid with '
               'keto end group(s) (hydroxylated, keto, epoxidated)\n'
               ' * SMILES: '
               'CC(\\C=C\\C1=C(C)CCCC1(C)C)=C/C=C/C(C)=C/C(O)C(O)\\C=C(C)\\C=C\\C=C(C)\\C=C\\C1=C(C)CCCC1(C)C '
               "NAME: 15,15'-dihydroxy-beta-carotene REASON: CORRECT Likely "
               'carotenoid based on carbon count and conjugation\n'
               ' * SMILES: '
               'C1(=C(C[C@H](CC1(C)C)O)C)/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C2=C(C[C@@H](CC2(C)C)O)C)\\C)\\C)/C)/C '
               'NAME: meso-zeaxanthin REASON: CORRECT Likely carotenoid based '
               'on carbon count and conjugation\n'
               ' * SMILES: '
               'O1C2(C1(CCCC2(C)C)C)/C=C/C(/C)=C/C=C/C(/C)=C/C=C/C=C(/C=C/C=C(/C=C/C=3C(CC(O)CC3C)(C)C)\\C)\\C '
               "NAME: 5',6'-Epoxy-5',6'-dihydro-beta,beta-caroten-3-ol REASON: "
               'CORRECT Likely carotenoid based on carbon count and '
               'conjugation\n'
               ' * SMILES: '
               'OC([C@@H](CC=C(C)C)CC/C(=C/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C1=C(CC[C@@H](C1(C)C)C/C=C(/CO)\\C)C)\\C)\\C)/C)/C)/C)(C)C '
               'NAME: '
               "(all-E,2R,2'R)-2-(4-hydroxy-3-methylbut-2-enyl)-2'-(3-methylbut-2-enyl)-1',2'-dihydro-beta,psi-caroten-1'-ol "
               'REASON: CORRECT Carotenoid with acyclic end group(s) '
               '(hydroxylated)\n'
               ' * SMILES: '
               'O=C1C(=C(/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C2=C(C(=O)[C@@H](O)CC2(C)C)C)\\C)\\C)/C)/C)C(C)(C)C[C@@H]1OC3OC(C(O)C(C3O)O)CO)C '
               "NAME: (3S,3'S)-Astaxanthin-beta-D-glucoside REASON: CORRECT "
               'Carotenoid with keto end group(s) (hydroxylated, keto)\n'
               ' * SMILES: '
               'CC(\\C=C\\C1=C(C)CCCC1(C)C)=C/C=C/C(C)=C/C1OC1\\C=C(C)\\C=C\\C=C(C)\\C=C\\C1=C(C)CCCC1(C)C '
               "NAME: beta-carotene 15,15'-epoxide REASON: CORRECT Likely "
               'carotenoid based on carbon count and conjugation\n'
               ' * SMILES: '
               'O[C@@]1(C(C[C@H](O)C[C@]1(O)C)(C)C)/C=C/C(/C)=C/C=C/C(/C)=C/C=C/C=C(/C=C/C=C(/C=C/[C@@]2(O)C(C[C@H](O)C[C@]2(O)C)(C)C)\\C)\\C '
               'NAME: Mactraxanthin REASON: CORRECT Likely carotenoid based on '
               'carbon count and conjugation\n'
               ' * SMILES: '
               'CC(\\C=C\\C=C(C)\\C=C\\C=C(C)\\C=C\\[C@H](O)C(C)(C)O)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/C1=C(C)CCCC1(C)C '
               'NAME: Deoxymyxol REASON: CORRECT Likely carotenoid based on '
               'carbon count and conjugation\n'
               ' * SMILES: '
               'O=C1C(=C(/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C=C(/C=C/C(O)C(O[C@@H]2OC[C@@H](O)[C@@H]([C@H]2O)O)(C)C)\\C)\\C)\\C)/C)/C)C(C)(C)CC1O)C '
               "NAME: 1'-beta-xylosyl-3,2'-dihydroxy-beta,psi-caroten-4-one "
               'REASON: CORRECT Carotenoid with keto end group(s) '
               '(hydroxylated, keto)\n'
               ' * SMILES: '
               'O[C@H]1CC([C@H](C(=C1)C)\\C=C\\C(\\C)=C\\C=C\\C(\\C)=C\\C=C\\C=C(\\C=C\\C=C(\\C=C\\C=2C(C[C@@H](O)CC2C)(C)C)/CO)/C)(C)C '
               'NAME: loroxanthin REASON: CORRECT Carotenoid with gamma end '
               'group(s) (hydroxylated)\n'
               ' * SMILES: '
               'OC/C(=C/CC1C(C(C(C)=CC1)/C=C/C(=C/C=C/C(=C/C=C/C=C(/CC/C=C(/CC/C=C(/CCC=C(C)C)\\C)\\C)\\C)/C)/C)(C)C)/C '
               'NAME: Nonaprenoxanthin REASON: CORRECT Carotenoid with gamma '
               'end group(s) (hydroxylated)\n'
               ' * SMILES: '
               'CC(\\C=C\\C=C(C)\\C=C\\C1=C(C)C[C@@H](O)CC1(C)C)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/C1=C(C)C[C@@H](O)CC1(C)C '
               'NAME: zeaxanthin REASON: CORRECT Likely carotenoid based on '
               'carbon count and conjugation\n'
               ' * SMILES: '
               'C1(=CCC(CC=C(C)C)C(C1/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C=C(/CCC=C(C)C)\\C)\\C)\\C)/C)/C)(C)C)C '
               'NAME: Nonaprene REASON: CORRECT Carotenoid with gamma end '
               'group(s)\n'
               ' * SMILES: '
               'CC(\\C=C\\C=C(C)\\C=C\\c1c(C)ccc(C)c1C)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/C1=C(C)CCCC1(C)C '
               'NAME: beta-Isorenieratene REASON: CORRECT Likely carotenoid '
               'based on carbon count and conjugation\n'
               ' * SMILES: '
               'OC1CC(C(C(=C1)C)\\C=C\\C(\\C)=C\\C=C\\C(\\C)=C\\C=C\\C=C(\\C=C\\C=C(\\C=C\\C2C(CC(=O)C=C2C)(C)C)/C)/C)(C)C '
               "NAME: 3'-hydroxy-e,e-caroten-3-one REASON: CORRECT Carotenoid "
               'with gamma, keto end group(s) (hydroxylated, keto)\n'
               ' * SMILES: '
               'O=C(O)/C(=C/C=C/C(=C/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C1=C(CC[C@H](C1(C)C)O)C)\\C)\\C)/C)/C)/C)/C '
               'NAME: 2-hydroxytorularhodin REASON: CORRECT Likely carotenoid '
               'based on carbon count and conjugation\n'
               ' * SMILES: '
               'O([C@H]1CC(C(=C(C1=O)C)/C=C/C(/C)=C/C=C/C(/C)=C/C=C/C=C(/C=C/C=C(/C=C/C=2C(C[C@H](O)C(=O)C2C)(C)C)\\C)\\C)(C)C)[C@@H]3OC([C@@H](O)C(O)C3O)CO '
               'NAME: Astaxanthin glucoside REASON: CORRECT Carotenoid with '
               'keto end group(s) (hydroxylated, keto)\n'
               ' * SMILES: '
               'O=C/C(=C/C=C/C(=C/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C1=C(CCCC1(C)C)C)\\C)\\C)/C)/C)/C)/C '
               'NAME: Torularhodinaldehyde REASON: CORRECT Likely carotenoid '
               'based on carbon count and conjugation\n'
               'False positives: SMILES: '
               'O=C(OC)C(=CC=CC(=CC=CC(=CC=CC=C(C=CC=C(C=CCC(CC1OC(C(O)C(C1O)O)CO)(C)C)C)C)C)C)C '
               'NAME: methyl 5-glucosyl-5,6-dihydro-apo-4,4′-lycopenoate '
               'REASON: WRONGLY CLASSIFIED Likely carotenoid based on carbon '
               'count and conjugation\n'
               ' * SMILES: '
               'C(\\C=C\\C=C(C)C)(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C=C(/C(O)=O)\\C)\\C)\\C)/C)/C '
               "NAME: 4,4'-diapolycopen-4-oic acid REASON: WRONGLY CLASSIFIED "
               'Likely carotenoid based on carbon count and conjugation\n'
               ' * SMILES: '
               'O=C1NC[C@@H](C=CC=CC=CC=CC=CC(C[C@@H](CC=C(C=CC=CC=C1)C)O)=O)C '
               'NAME: Bombyxamycin A REASON: WRONGLY CLASSIFIED Carotenoid '
               'with keto end group(s) (hydroxylated, keto)\n'
               ' * SMILES: '
               'O=C(OC/C(=C/C=C/C(=C/C=C/C(=C/C=C/C(=C/C=C/C=C(/CCCCCO)\\C)/C)/C)/C)/C)C=1C2=C(C=CC=C2)NC1CO '
               'NAME: Rhodestrin REASON: WRONGLY CLASSIFIED Likely carotenoid '
               'based on carbon count and conjugation\n'
               ' * SMILES: '
               'CC(C)=CCC\\C(C)=C\\C=C\\C(\\C)=C\\C=C\\C(\\C)=C\\C=C/C=C(\\C)/C=C/C=C(\\C)/C=C/C=C(\\C)CCC=C(C)C '
               'NAME: 15-cis-lycopene REASON: WRONGLY CLASSIFIED Likely '
               'carotenoid based on carbon count and conjugation\n'
               ' * SMILES: '
               'CC(C)=CCC\\C(C)=C/C=C/C(/C)=C/C=C/C(/C)=C/C=C/C=C(\\C)/C=C/C=C(\\C)/C=C/C=C(\\C)CCC=C(C)C '
               'NAME: 5-cis-lycopene REASON: WRONGLY CLASSIFIED Likely '
               'carotenoid based on carbon count and conjugation\n'
               ' * SMILES: '
               'CC(C)=CCC\\C(C)=C\\CC\\C(C)=C/C=C/C(C)=C/C=C/C=C(C)/C=C/C=C(/C)CC\\C=C(/C)CCC=C(C)C '
               "NAME: 9,9'-di-cis-zeta-carotene REASON: WRONGLY CLASSIFIED "
               'Likely carotenoid based on carbon count and conjugation\n'
               ' * SMILES: '
               'O=C(O)C\\C=C\\C=C/C=C/C=C/C=C/C=C/C=C\\C=C\\C=C(\\C(=O)C)/C '
               'NAME: '
               '(3E,5Z,7E,9E,11E,13E,15Z,17E,19E)-20-methyl-21-oxodocosa-3,5,7,9,11,13,15,17,19-nonaenoic '
               'acid REASON: WRONGLY CLASSIFIED Carotenoid with keto end '
               'group(s) (hydroxylated, keto, carboxylated)\n'
               ' * SMILES: '
               'O(C(=O)/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(\\C=C\\C=1C(CCCC1C)(C)C)/C)\\C)/C)/C)CC '
               'NAME: all-trans-Carophyll yellow REASON: WRONGLY CLASSIFIED '
               'Likely carotenoid based on carbon count and conjugation\n'
               ' * SMILES: '
               'CC(C)=CCC\\C(C)=C\\C=C\\C(C)=C\\C=C\\C(C)=C\\C=C\\C=C(C)\\C=C\\C=C(C)\\C=C\\C=C(/C)CCC=C(C)C '
               'NAME: lycopene REASON: WRONGLY CLASSIFIED Likely carotenoid '
               'based on carbon count and conjugation\n'
               ' * SMILES: '
               'O=C1N[C@H](CC=C(C=CC=C(C=CC=C[C@@H]([C@@](C=CC=CC=C(C=C1OC)C)(O)C)O[C@@H]2OC[C@@H](O[C@@H]3O[C@@H]([C@@H](NC)CC3)C)[C@@H]([C@H]2NC)O)C)C)C '
               'NAME: Incednine REASON: WRONGLY CLASSIFIED Likely carotenoid '
               'based on carbon count and conjugation\n'
               ' * SMILES: '
               'C(/C=C/C=C(C)C)(=C\\C=C\\C(=C\\C=C\\C=C(\\C=C\\C=C(\\C=C\\C=C(\\C([O-])=O)/C)/C)/C)\\C)\\C '
               "NAME: 4,4'-diapolycopen-4-oate REASON: WRONGLY CLASSIFIED "
               'Likely carotenoid based on carbon count and conjugation\n'
               ' * SMILES: '
               'O=C(O[C@@H]1OC[C@@H](O)[C@@H]([C@H]1OC(=O)CCCCCCCCCCC(C)C)O)/C(=C/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C=C(/C(=O)O)\\C)\\C)\\C)/C)/C)/C '
               'NAME: Diapolycopenedioic acid xylosyl ester REASON: WRONGLY '
               'CLASSIFIED Likely carotenoid based on carbon count and '
               'conjugation\n'
               ' * SMILES: '
               'C\\C(CCCC(C)(C)O)=C/C=C/C(C)=C/C=C/C(C)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/c1c(C)ccc(C)c1C '
               'NAME: Hydroxychlorobactene REASON: WRONGLY CLASSIFIED Likely '
               'carotenoid based on carbon count and conjugation\n'
               ' * SMILES: '
               'O[C@H]1CC(C(=C(C1=O)C)/C=C/C(/C)=C/C=C/C(/C)=C/C=C/C=C(/C=C/C=C(/C#CC=2C(C[C@H](O)C(=O)C2C)(C)C)\\C)\\C)(C)C '
               'NAME: 7,8-Didehydroastaxanthin REASON: WRONGLY CLASSIFIED '
               'Carotenoid with keto end group(s) (hydroxylated, keto)\n'
               ' * SMILES: '
               'O[C@H]1CC(\\C(\\C(=C1)C)=C\\C=C(/C)\\C=C\\C=C(/C)\\C=C\\C=C\\C(=C\\C=C\\C(=C\\C=C/2\\C(CC(=O)C=C2C)(C)C)\\C)\\C)(C)C '
               'NAME: eschscholtzxanthone REASON: WRONGLY CLASSIFIED '
               'Carotenoid with gamma, keto end group(s) (hydroxylated, keto)\n'
               ' * SMILES: '
               'C\\C(CO)=C/C=C/C(C)=C/C=C/C(C)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/C1=C(C)CCCC1(C)C '
               "NAME: 4'-apo-beta-carotenol REASON: WRONGLY CLASSIFIED Likely "
               'carotenoid based on carbon count and conjugation\n'
               ' * SMILES: '
               'CC(C)=C\\C=C\\C(C)=C\\C=C\\C(C)=C\\C=C\\C(C)=C\\C=C\\C=C(C)\\C=C\\C=C(C)\\C=C\\C1=C(C)CCCC1(C)C '
               'NAME: torulene REASON: WRONGLY CLASSIFIED Likely carotenoid '
               'based on carbon count and conjugation\n'
               ' * SMILES: '
               'O=C(OC)/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C=C(/C=C/CC(OC1OC(C(O)C(C1O)O)CO)(C)C)\\C)\\C)\\C)/C)/C '
               'NAME: methyl glucosyl-3,4-dehydro-apo-8′-lycopenoate REASON: '
               'WRONGLY CLASSIFIED Likely carotenoid based on carbon count and '
               'conjugation\n'
               ' * SMILES: '
               'O[C@H]1CC(C(=C(C1)C)C#C/C(/C)=C/C=C/C(/C)=C/C=C/C=C(/C=C/C=C(/C#CC=2C(C[C@H](O)CC2C)(C)C)\\C)\\C)(C)C '
               'NAME: '
               'Alloxanthin/Tetradehydrozeaxanthin/(Cynthiaxanthin)/(Pectenoxanthin) '
               'REASON: WRONGLY CLASSIFIED Likely carotenoid based on carbon '
               'count and conjugation\n'
               ' * SMILES: '
               'CC(C)=CCC\\C(C)=C\\C=C\\C(C)=C\\C=C\\C(C)=C\\C=C\\C=C(C)\\C=C\\C=C(C)\\C=C\\C=C(C)\\C=C\\C=C(C)C '
               'NAME: all-trans-3,4-didehydrolycopene REASON: WRONGLY '
               'CLASSIFIED Likely carotenoid based on carbon count and '
               'conjugation\n'
               ' * SMILES: '
               'ClC1=C2C(C(=O)[C@](C1=O)(OC(=O)C)C)=CN(CCCCN3C(=CC4=C(Cl)C(=O)[C@@](C(C4=C3)=O)(OC(=O)C)C)C=CC(=C[C@H](CC)C)C)C(=C2)C=CC(=C[C@H](CC)C)C '
               'NAME: Sclerotiorin E REASON: WRONGLY CLASSIFIED Carotenoid '
               'with keto end group(s) (hydroxylated, keto, carboxylated)\n'
               ' * SMILES: '
               'O=C(O)/C=C/C/C=C/C=C/C=C\\C=C\\C=C/C=C/C=C/C=C\\C=C\\C=C(\\C(=O)C)/C '
               'NAME: 2-Dehydro-3-deoxylaetiporic acid A REASON: WRONGLY '
               'CLASSIFIED Carotenoid with keto end group(s) (hydroxylated, '
               'keto, carboxylated)\n'
               ' * SMILES: '
               'O1[C@@H](OC(CCC/C(/C)=C/C=C/C(/C)=C/C=C/C(/C)=C/C=C/C=C(/C=C/C=C(/C=C/C=2C(=C(C=CC2C)C)C)\\C)\\C)(C)C)C(O)C(O)[C@H](O)[C@H]1CO '
               'NAME: OH-Chlorobactene glucoside REASON: WRONGLY CLASSIFIED '
               'Likely carotenoid based on carbon count and conjugation\n'
               ' * SMILES: '
               'O[C@H]1CC(C(=C(C1)C)\\C=C\\C(\\C)=C\\C=C\\C(\\C)=C\\C=C\\C=C(\\C=C\\C=C(\\C#CC=2C(C[C@H](O)CC2C)(C)C)/C)/C)(C)C '
               'NAME: Diatoxanthin REASON: WRONGLY CLASSIFIED Likely '
               'carotenoid based on carbon count and conjugation\n'
               'False negatives: SMILES: '
               'CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\[C@H]1[C@H](COP(O)(=O)OP(O)(O)=O)[C@@]1(C)CC\\C=C(/C)CC\\C=C(/C)CCC=C(C)C '
               'NAME: (1S,2S,3S)-prephytoene diphosphate REASON: MISSED No '
               'characteristic polyene chain found\n'
               ' * SMILES: '
               'O=C(OCC1OC(O[C@@H]2CC(=C(/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C3=C(C[C@@H](OC4C(O)C(O)C(O)C(C4)COC(=O)CCCCCCCCCC(C)C)CC3(C)C)C)\\C)\\C)/C)/C)C(C2)(C)C)C)C(O)C(C1O)O)CCCCCCCCCC(C)C '
               'NAME: Thermobiszeaxanthin Z2-13-13 REASON: MISSED Too many '
               'carbons (79) for a carotenoid\n'
               ' * SMILES: '
               '[H]C(=O)\\C=C\\C=C(C)\\C=C\\C=C(C)\\C=C\\C1=C(C)CCCC1(C)C '
               "NAME: 14'-apo-beta-carotenal REASON: MISSED Insufficient "
               'conjugation (6 double bonds)\n'
               ' * SMILES: '
               'CC(\\C=C\\C=C(/C)C(O)=O)=C/C=C/C=C(C)/C=C/C=C(\\C)C(O)=O NAME: '
               'crocetin REASON: MISSED No characteristic end groups found\n'
               ' * SMILES: '
               'C(\\C=C\\C(OC)=O)(=C\\C=C\\C(=C\\C=C\\C=C(\\C=C\\C=C(\\C=C\\C(=O)O)/C)/C)\\C)/C '
               'NAME: bixin REASON: MISSED No characteristic end groups found\n'
               ' * SMILES: '
               'O=C(OC[C@H]1O[C@@H](OC(C/C=C/C(=C/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C=C(/C=C/CC(O)(C)C)\\C)\\C)\\C)/C)/C)/C)(C)C)[C@H](O)[C@H]([C@@H]1O)O)CCCCCCCCCCCCC '
               'NAME: I\' -(6"-0-myristyl-p-glucopyanosyloxy) -3,4,3\', 4\' '
               "-tetradehydro1,2,1',2'-tetrahydro-l-hydroxy-yl,yl-carotene "
               'REASON: MISSED Too many carbons (60) for a carotenoid\n'
               ' * SMILES: '
               'OC1(C(C[C@@H](O)C=C1C)(C)C)/C=C/C(/C)=C/C=C/C(/C)=C/C=C/C=C(/C=C/C=C(/C(=O)CC=2C(C[C@H](O)CC2C)(C)C)\\COC(=O)CCCCCCCCCCCC)\\C '
               "NAME: 6'-Hydroxysiphonaxanthin ester/ 6'-Hydroxysiphonaxanthin "
               "dodecenoate/ 6'-Hydroxysiphonaxanthin "
               '19-(trans-delta2-dodecenoate) REASON: MISSED Too many carbons '
               '(53) for a carotenoid\n'
               ' * SMILES: '
               'CC(\\C=C\\C=C(C)\\C=C\\C1=C(C)C[C@H](CC1(C)C)O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/C1=C(C)C[C@H](CC1(C)C)O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O '
               'NAME: zeaxanthin bis(beta-D-glucoside) REASON: MISSED Too many '
               'carbons (52) for a carotenoid\n'
               ' * SMILES: '
               'O[C@@H]1CC([C@H](C(=C1)C)\\C=C\\C(\\C)=C\\C=C\\C(\\C)=C\\C=C\\C=C(\\C=C\\C=C(\\C(=O)CC=2C(C[C@H](O)CC2C)(C)C)/COC(=O)/C=C/CCCCCCCCC)/C)(C)C '
               'NAME: Siphonaxanthin ester/ Siphonaxanthin dodecenoate/ '
               '(Siphonein) REASON: MISSED Too many carbons (52) for a '
               'carotenoid\n'
               ' * SMILES: '
               'CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\[C@@H]1[C@@H](COP(O)(=O)OP(O)(O)=O)[C@]1(C)CC\\C=C(/C)CC\\C=C(/C)CCC=C(C)C '
               'NAME: (1R,2R,3R)-prephytoene diphosphate REASON: MISSED No '
               'characteristic polyene chain found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C[C@H]1CN(C(=O)C2=C(C3=CC=CC=C3CO[C@@H]1CN(C)C(=O)NC4=CC=C(C=C4)OC)C5=CC=CC=C5N2C)[C@H](C)CO',
                                     'name': 'LSM-30659',
                                     'reason': 'No characteristic carotenoid '
                                               'backbone found'},
                                 {   'smiles': 'c1c[nH+]c[nH]1',
                                     'name': 'imidazolium cation',
                                     'reason': 'Too few carbons (3) for a '
                                               'carotenoid'},
                                 {   'smiles': 'ClC=1C(O)=C2O[C@@]([C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC)[C@@H](O)C(=C2)C1)C(C)=C)C(=O)N3[C@@H](C=CC3)C(=O)N/C(=C(/CC)\\C)/C(=O)N/C(=C/C(O)=O)/C(O)=O)(CC)C',
                                     'name': 'Phomopsin A',
                                     'reason': 'No characteristic carotenoid '
                                               'backbone found'},
                                 {   'smiles': 'ClC(Cl)[C@H](O)CC=1OC(=O)C=2C(O)=CC(=CC2C1)O',
                                     'name': 'Desmethyldichlorodiaportin',
                                     'reason': 'Too few carbons (12) for a '
                                               'carotenoid'},
                                 {   'smiles': 'ClC=1C(=C(O)C2=C(C1O)C(=O)C=CC2=O)CC=C(C)C',
                                     'name': 'Chlorosesamone',
                                     'reason': 'Too few carbons (15) for a '
                                               'carotenoid'},
                                 {   'smiles': 'CCCCC[C@H]1O[C@H]1C\\C=C/CCCCCCCC(O)=O',
                                     'name': '(+)-vernolic acid',
                                     'reason': 'Too few carbons (18) for a '
                                               'carotenoid'},
                                 {   'smiles': 'O(C1C(O)C(OC(OC=2C3=C(C=CC2)C=C(C(=C3O)C(=O)C)C)C1O)CO)C4OC(C(O)C(O)C4O)CO',
                                     'name': 'Orientaloside',
                                     'reason': 'No characteristic carotenoid '
                                               'backbone found'},
                                 {   'smiles': 'O(C1=C(C=CC=C1C(O)=O)C)C(=O)C',
                                     'name': 'CRESOPYRINE',
                                     'reason': 'Too few carbons (10) for a '
                                               'carotenoid'},
                                 {   'smiles': 'C[C@H]1CN(S(=O)(=O)C2=C(C=C(C=C2)C#CC3=CC=NC=C3)O[C@@H]1CN(C)C(=O)NC4=CC5=C(C=C4)OCO5)[C@@H](C)CO',
                                     'name': '3-(1,3-benzodioxol-5-yl)-1-[[(4S,5S)-2-[(2S)-1-hydroxypropan-2-yl]-4-methyl-1,1-dioxo-8-(2-pyridin-4-ylethynyl)-4,5-dihydro-3H-6,1$l^{6},2-benzoxathiazocin-5-yl]methyl]-1-methylurea',
                                     'reason': 'No characteristic carotenoid '
                                               'backbone found'},
                                 {   'smiles': 'O1C2C3C(CCC3=C)C(CCC2C(C1=O)=C)=C',
                                     'name': '3,6,9-Trimethylidene-3a,4,5,6a,7,8,9a,9b-octahydroazuleno[4,5-b]furan-2-one',
                                     'reason': 'Too few carbons (15) for a '
                                               'carotenoid'}],
    'sample_false_negatives': [   {   'smiles': 'CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\[C@H]1[C@H](COP(O)(=O)OP(O)(O)=O)[C@@]1(C)CC\\C=C(/C)CC\\C=C(/C)CCC=C(C)C',
                                      'name': '(1S,2S,3S)-prephytoene '
                                              'diphosphate',
                                      'reason': 'Does not match carotenoid '
                                                'structural requirements'},
                                  {   'smiles': 'O=C(OCC1OC(O[C@@H]2CC(=C(/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C3=C(C[C@@H](OC4C(O)C(O)C(O)C(C4)COC(=O)CCCCCCCCCC(C)C)CC3(C)C)C)\\C)\\C)/C)/C)C(C2)(C)C)C)C(O)C(C1O)O)CCCCCCCCCC(C)C',
                                      'name': 'Thermobiszeaxanthin Z2-13-13',
                                      'reason': 'Too many carbons (79) for a '
                                                'carotenoid'},
                                  {   'smiles': 'CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\[C@@H]1[C@@H](COP(O)(=O)OP(O)(O)=O)[C@]1(C)CC\\C=C(/C)CC\\C=C(/C)CCC=C(C)C',
                                      'name': '(1R,2R,3R)-prephytoene '
                                              'diphosphate',
                                      'reason': 'Does not match carotenoid '
                                                'structural requirements'},
                                  {   'smiles': 'O=C1C(=C(/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C=C(/CCCC(OC2OC(C(O)C(C2O)O)COC(=O)C(C(O)CCCCCCCCCCCCCCCCCCCC(/C=C\\CCCCCCCCCCCCCC(=O)O)C)CCCCCCCCCCCCCCCCCCCC)(C)C)\\C)\\C)\\C)/C)/C)C(C)(C)CC1)C',
                                      'name': "1'-[(6-O-Mycoloyl-beta-D-glucopyranosyl)",
                                      'reason': 'Too many carbons (106) for a '
                                                'carotenoid'},
                                  {   'smiles': 'CCCCCCCCCCCCCCCC(=O)O[C@@H]1CC(C)=C(\\C=C\\C(C)=C\\C=C\\C(C)=C\\C=C\\C=C(C)\\C=C\\C=C(C)\\C=C\\[C@H]2C(C)=C[C@@H](CC2(C)C)OC(=O)CCCCCCCCCCCCCCC)C(C)(C)C1',
                                      'name': 'Helenien',
                                      'reason': 'Too many carbons (72) for a '
                                                'carotenoid'},
                                  {   'smiles': 'CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\C1C(COP(O)(=O)OP(O)(O)=O)C1(C)CC\\C=C(/C)CC\\C=C(/C)CCC=C(C)C',
                                      'name': 'prephytoene diphosphate',
                                      'reason': 'Does not match carotenoid '
                                                'structural requirements'},
                                  {   'smiles': 'O=C1C(=C(/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C=C(/C=C/[C@H](O)C(O[C@H]2O[C@H]([C@H](O)[C@H]([C@@H]2O)O)COC(=O)CCCCCCCCCCCC(C)C)(C)C)\\C)\\C)\\C)/C)/C)C(C)(C)CC1)C',
                                      'name': 'Salinixanthin',
                                      'reason': 'Too many carbons (61) for a '
                                                'carotenoid'},
                                  {   'smiles': 'C(C\\C=C(\\CCC=C(C)C)/C)/C(=C/CC/C(/C)=C\\C=C\\C=C(/C=C\\C=C(\\CC\\C=C(\\CCC=C(C)C)/C)/C)\\C)/C',
                                      'name': "(15Z,9'Z)-7,7',8,8',11,12-Hexahydrolycopene",
                                      'reason': 'Does not match carotenoid '
                                                'structural requirements'},
                                  {   'smiles': 'O=C(OCC1OC(O[C@@H]2CC(=C(/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C3=C(CCCC3(C)C)C)\\C)\\C)/C)/C)C(C2)(C)C)C)C(O)C(C1O)O)CCCCCCCCCCCC(C)C',
                                      'name': 'C1-15 thermocryptoxanthin-15',
                                      'reason': 'Too many carbons (61) for a '
                                                'carotenoid'},
                                  {   'smiles': 'O=C1C(=C(/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C2=C(C(=O)[C@@H](OC3OC(C(O)C(C3O)O)C)CC2(C)C)C)\\C)\\C)/C)/C)C(C)(C)C[C@@H]1OC4OC(C(O)C(C4O)O)C)C',
                                      'name': 'Astaxanthin dirhamnoside',
                                      'reason': 'Too many carbons (52) for a '
                                                'carotenoid'}],
    'attempt': 2,
    'reasoning': 'Let me analyze the issues from the previous attempt:\n'
                 '\n'
                 '1. False Positives Analysis:\n'
                 '- Many molecules with extended conjugation but not true '
                 'carotenoids were misclassified\n'
                 '- The carbon count range (20-50) was too permissive\n'
                 '- The program relied too heavily on conjugation patterns '
                 'without considering the full structure\n'
                 '\n'
                 '2. False Negatives Analysis:\n'
                 '- Some legitimate carotenoids were missed due to strict '
                 'carbon count limits (>50 C)\n'
                 '- The polyene chain detection was too rigid\n'
                 '- Some carotenoid precursors (like prephytoene diphosphate) '
                 'were missed\n'
                 '- Modified carotenoids with sugar groups were excluded\n'
                 '\n'
                 "Here's how we can improve:\n"
                 '\n'
                 '1. Better Carbon Count Logic:\n'
                 '- Keep base range around 40±10 for typical carotenoids\n'
                 '- Add special exceptions for known modifications '
                 '(glycosides, esters)\n'
                 '- Consider molecular weight as additional filter\n'
                 '\n'
                 '2. Improved Structural Detection:\n'
                 '- Add more specific patterns for carotenoid backbone\n'
                 '- Better handle branching methyl groups\n'
                 '- Include patterns for carotenoid precursors\n'
                 '- Add more end group patterns\n'
                 '\n'
                 '3. Additional Checks:\n'
                 '- Count methyl branches\n'
                 '- Verify spacing of conjugated system\n'
                 '- Check for characteristic carbon skeleton\n'
                 '- Add molecular symmetry check\n'
                 '\n'
                 "Here's the improved code:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 190,
    'num_false_positives': 100,
    'num_true_negatives': 40736,
    'num_false_negatives': 11,
    'num_negatives': None,
    'precision': 0.6551724137931034,
    'recall': 0.945273631840796,
    'f1': 0.7739307535641547,
    'accuracy': 0.9972951239125667,
    'negative_predictive_value': 0.9997300414754461}