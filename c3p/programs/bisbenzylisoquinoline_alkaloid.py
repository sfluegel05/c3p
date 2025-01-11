"""
Classifies: CHEBI:133004 bisbenzylisoquinoline alkaloid
"""
"""
Classifies: CHEBI:bisbenzylisoquinoline alkaloid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_bisbenzylisoquinoline_alkaloid(smiles: str):
    """
    Determines if a molecule is a bisbenzylisoquinoline alkaloid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a bisbenzylisoquinoline alkaloid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for minimum complexity
    if mol.GetNumAtoms() < 30:  # These are large molecules
        return False, "Molecule too small for bisbenzylisoquinoline alkaloid"

    # Look for isoquinoline cores (can be N-methylated or not)
    isoquinoline_pattern = Chem.MolFromSmarts("[#6]1[#6][#6][#7][#6][#6]2[#6](=[#6][#6]=1)[#6]=[#6][#6]=2")
    isoquinoline_matches = mol.GetSubstructMatches(isoquinoline_pattern)
    
    if len(isoquinoline_matches) < 2:
        return False, "Must contain two isoquinoline cores"

    # Look for ether bridges (-O-)
    ether_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6]")
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    
    if len(ether_matches) < 1:
        return False, "Must contain ether bridges between units"

    # Check for N-methyl groups or NH (typical in these alkaloids)
    n_methyl_pattern = Chem.MolFromSmarts("[#7]-[CH3]")
    nh_pattern = Chem.MolFromSmarts("[#7;H1]")
    n_methyl_matches = mol.GetSubstructMatches(n_methyl_pattern)
    nh_matches = mol.GetSubstructMatches(nh_pattern)
    
    if len(n_methyl_matches) + len(nh_matches) < 2:
        return False, "Must contain two nitrogen centers (N-methyl or NH)"

    # Look for benzyl groups attached to isoquinoline
    benzyl_pattern = Chem.MolFromSmarts("[#6]1[#6][#6][#6][#6][#6]1-[CH2]")
    benzyl_matches = mol.GetSubstructMatches(benzyl_pattern)
    
    if len(benzyl_matches) < 2:
        return False, "Must contain two benzyl groups"

    # Optional: Check for common substituents (methoxy, hydroxy)
    methoxy_pattern = Chem.MolFromSmarts("[#6]-[#8]-[CH3]")
    hydroxy_pattern = Chem.MolFromSmarts("[#6]-[#8;H1]")
    methoxy_matches = mol.GetSubstructMatches(methoxy_pattern)
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    
    # Count aromatic rings
    aromatic_rings = 0
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            aromatic_rings += 1
            
    if aromatic_rings < 4:  # Should have at least 4 aromatic rings
        return False, "Insufficient number of aromatic rings"

    # Additional check for methylenedioxy groups (common in this class)
    methylenedioxy_pattern = Chem.MolFromSmarts("[#6]1[#8][#6][#8]1")
    methylenedioxy_matches = mol.GetSubstructMatches(methylenedioxy_pattern)

    # Verify overall connectivity
    if not mol.HasSubstructMatch(Chem.MolFromSmarts("[#6]1~[#6]~[#6]~[#7]~[#6]~[#6]~[#6]~[#6]~1")):
        return False, "Missing key ring connectivity"

    return True, "Contains two benzylisoquinoline units connected by ether bridges with appropriate substitution patterns"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:133004',
                          'name': 'bisbenzylisoquinoline alkaloid',
                          'definition': 'A type of benzylisoquinoline alkaloid '
                                        'whose structures are built up of two '
                                        'benzylisoquinoline units linked by '
                                        'ether bridges. Various structural '
                                        'patterns resulting from additional '
                                        'bridging between the two units by '
                                        'direct carbon-carbon bridging or by '
                                        'methylenedioxy groups are common.',
                          'parents': ['CHEBI:22750'],
                          'xrefs': [   'PMID:1955879',
                                       'PMID:2191354',
                                       'PMID:3323421'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'N12CN3CN(CN(C1)CC3)CC2',
                                     'name': '1,3,6,8-tetraazatricyclo[4,4,1,1(3,8)]dodecane',
                                     'reason': 'Molecule too small for '
                                               'bisbenzylisoquinoline '
                                               'alkaloid'},
                                 {   'smiles': 'CC1CC2C(CC=C1CCO)C(C)C(=O)O2',
                                     'name': '(8alpha,10beta,11beta)-3-hydroxy-4,15-dinor-1(5)-xanthen-12,8-olide',
                                     'reason': 'Molecule too small for '
                                               'bisbenzylisoquinoline '
                                               'alkaloid'},
                                 {   'smiles': 'CCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](COP(O)(O)=O)[C@H](O)\\C=C\\CCCCCCCCCCCCC',
                                     'name': 'N-docosanoylsphingosine-1-phosphate',
                                     'reason': 'Must contain two isoquinoline '
                                               'cores'},
                                 {   'smiles': 'C1CNCCC1SC2=NC(=CC=C2)Cl',
                                     'name': '2-chloro-6-(4-piperidinylthio)pyridine',
                                     'reason': 'Molecule too small for '
                                               'bisbenzylisoquinoline '
                                               'alkaloid'},
                                 {   'smiles': '[H][C@@]1(CC[C@@]2(C)[C@]3([H])CC[C@@]4([H])[C@@](C)([C@H](CC[C@@]44C[C@@]34CC[C@]12C)OC(C)=O)C(O)=O)[C@H](C)CCC(=O)C(C)C',
                                     'name': 'bonianic acid B, (rel)-',
                                     'reason': 'Must contain two isoquinoline '
                                               'cores'},
                                 {   'smiles': 'O=C1OC[C@@]2(C1=C(CC[C@]34[C@H]2C[C@H](CC[C@H]3C)C4(C)C)C)C',
                                     'name': 'Harzianelactone',
                                     'reason': 'Molecule too small for '
                                               'bisbenzylisoquinoline '
                                               'alkaloid'},
                                 {   'smiles': 'O=C(NC1=NC=C(Cl)C=C1)CSC2=NN=C3C(NC4=C3C=CC=C4)=N2',
                                     'name': 'dCeMM2',
                                     'reason': 'Molecule too small for '
                                               'bisbenzylisoquinoline '
                                               'alkaloid'},
                                 {   'smiles': 'N[C@@H](CS[C@H](\\C=C\\C=C\\C=C/C=C/CC(O)=O)[C@@H](O)CCCC(O)=O)C(O)=O',
                                     'name': '(13E)-16-carboxy-Delta(13)-17,18,19,20-tetranor-leukotriene '
                                             'E4',
                                     'reason': 'Molecule too small for '
                                               'bisbenzylisoquinoline '
                                               'alkaloid'},
                                 {   'smiles': 'C[C@@H]1O[C@@H](OC[C@H]2OC(O)[C@H](NC(C)=O)[C@@H](O[C@@H]3O[C@@H](C)[C@@H](O)[C@@H](O)[C@@H]3O)[C@@H]2O[C@@H]2O[C@H](CO)[C@@H](O[C@@H]3O[C@H](CO[C@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]4O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4NC(C)=O)[C@@H](O)[C@H](O[C@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]4O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4NC(C)=O)[C@@H]3O)[C@H](O)[C@H]2NC(C)=O)[C@@H](O)[C@H](O)[C@@H]1O',
                                     'name': 'alpha-L-Fucp-(1->3)-[alpha-L-Fucp-(1->6)]-{beta-D-GlcpNAc-(1->2)-alpha-D-Manp-(1->3)-[beta-D-GlcpNAc-(1->2)-alpha-D-Manp-(1->6)]-beta-D-Manp-(1->4)-beta-GlcpNAc-(1->4)}-D-GlcpNAc',
                                     'reason': 'Must contain two isoquinoline '
                                               'cores'},
                                 {   'smiles': 'CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1ccc(N)nc1=O)OC(=O)CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CC',
                                     'name': 'CDP-1-stearoyl-2-(4Z,7Z,10Z,13Z,16Z,19Z)-docosahexaenoyl-sn-glycerol',
                                     'reason': 'Must contain two isoquinoline '
                                               'cores'}],
    'sample_false_negatives': [   {   'smiles': 'COc1cc2CCN(C)[C@H](Cc3ccc(Oc4cc(C[C@H]5N(C)CCc6cc(OC)c(O)cc56)ccc4O)cc3)c2cc1O',
                                      'name': 'guattegaumerine',
                                      'reason': 'Must contain two isoquinoline '
                                                'cores'},
                                  {   'smiles': 'O1C=2C=C3C(=CC2OC)CCN([C@H]3CC4=CC=C(OC=5C=C(C[C@@]6(N(CCC=7C6=C1C(O)=C(OC)C7)C)[H])C=CC5OC)C=C4)C',
                                      'name': 'fangchinoline',
                                      'reason': 'Must contain two isoquinoline '
                                                'cores'},
                                  {   'smiles': 'COc1ccc2C[C@H]3[C@H]4C[C@]5([C@H]6Oc1c2[C@@]46CCN3C)N1CCc2ccc(OC)c3Oc4c(O)c(OC)ccc4C(C5=O)=C1c23',
                                      'name': 'Cancentrine',
                                      'reason': 'Must contain two isoquinoline '
                                                'cores'},
                                  {   'smiles': 'COc1cc2CCN=C3Cc4ccc(Oc5c(OC)ccc(C[C@H]6N(C)CCc7cc(OC)c(OC)c(Oc1cc23)c67)c5OC)cc4',
                                      'name': 'Calafatimine',
                                      'reason': 'Must contain two isoquinoline '
                                                'cores'},
                                  {   'smiles': 'COc1cc2CCN(C)[C@@H](Cc3ccc(Oc4cc(C[C@H]5N(C)CCc6cc(OC)c(O)cc56)ccc4O)cc3)c2cc1O',
                                      'name': 'berbamunine',
                                      'reason': 'Must contain two isoquinoline '
                                                'cores'},
                                  {   'smiles': 'COc1ccc2C[C@H]3N(C)CCc4cc(OC)c(Oc5c(OC)c(OC)cc6CCN(C)[C@@H](Cc7ccc(Oc1c2)cc7)c56)cc34',
                                      'name': 'Obaberine',
                                      'reason': 'Must contain two isoquinoline '
                                                'cores'},
                                  {   'smiles': 'COc1cc2CCN(C)[C@H](Cc3ccc(Oc4cc(C[C@H]5N(C)CCc6cc(OC)c(OC)cc56)ccc4O)cc3)c2cc1OC',
                                      'name': 'dauricine',
                                      'reason': 'Must contain two isoquinoline '
                                                'cores'},
                                  {   'smiles': 'COc1cc2CCN(C)[C@@H]3Cc4ccc(Oc5c(O)c(OC)cc6CCN(C)[C@H](Cc7ccc(Oc(c1O)c23)cc7)c56)cc4',
                                      'name': 'Isochondrodendrine',
                                      'reason': 'Must contain two isoquinoline '
                                                'cores'},
                                  {   'smiles': '[H][C@@]12Cc3ccc(O)c(Oc4ccc(C[C@]5([H])N(C)CCc6cc(OC)c(OC)c(Oc7cc1c(CCN2C)cc7OC)c56)cc4)c3',
                                      'name': 'oxyacanthine',
                                      'reason': 'Must contain two isoquinoline '
                                                'cores'},
                                  {   'smiles': 'COc1cc2CCN(C)[C@H]3Cc4ccc(O)c(Oc5ccc(C[C@H]6N(C)CCc7c6cc(Oc1cc23)c(O)c7OC)cc5)c4',
                                      'name': 'Aromoline',
                                      'reason': 'Must contain two isoquinoline '
                                                'cores'}],
    'attempt': 0,
    'reasoning': 'Let me help create a program to classify '
                 'bisbenzylisoquinoline alkaloids. These molecules have some '
                 'key structural features:\n'
                 '\n'
                 '1. Two benzylisoquinoline units\n'
                 '2. Units are connected by ether bridges (O-linkages)\n'
                 '3. May have additional C-C bridges or methylenedioxy groups\n'
                 '4. Typically contain N-methylated or N-unsubstituted '
                 'isoquinoline rings\n'
                 '5. Often have methoxy (OMe) and hydroxy (OH) substituents\n'
                 '\n'
                 "Here's the program:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 142268,
    'num_false_negatives': 32,
    'num_negatives': None,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': 0.9997751229796206,
    'negative_predictive_value': 0.9997751229796206}