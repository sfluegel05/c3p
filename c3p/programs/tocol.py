"""
Classifies: CHEBI:39437 tocol
"""
"""
Classifies: tocol compounds
Definition: A chromanol with a chroman-6-ol skeleton substituted at position 2 
by a saturated or triply-unsaturated hydrocarbon chain of three isoprenoid units
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tocol(smiles: str):
    """
    Determines if a molecule is a tocol based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (is_tocol, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for basic chroman core (benzene fused to pyran)
    chroman_pattern = Chem.MolFromSmarts("O1CCCc2ccccc12")
    if not mol.HasSubstructMatch(chroman_pattern):
        return False, "No chroman core found"
    
    # Check for hydroxyl group (or modified hydroxyl) at position 6
    # Note: Position 6 corresponds to para position relative to oxygen
    hydroxy_pattern = Chem.MolFromSmarts("O1CCCc2cc([OH,O])ccc12")
    hydroxy_ester_pattern = Chem.MolFromSmarts("O1CCCc2cc(O[C,P])ccc12")
    if not (mol.HasSubstructMatch(hydroxy_pattern) or mol.HasSubstructMatch(hydroxy_ester_pattern)):
        return False, "No hydroxyl group (or derivative) at position 6"
    
    # Check for substitution at position 2 (carbon next to oxygen in pyran ring)
    # Looking for branched carbon chain
    c2_substitution = Chem.MolFromSmarts("O1[C](CC)CCc2ccccc12")
    if not mol.HasSubstructMatch(c2_substitution):
        return False, "No substitution at position 2"
    
    # Count carbons in the longest chain from position 2
    # This checks for the presence of the isoprenoid chain
    side_chain_pattern = Chem.MolFromSmarts("[CH3]C(C)CCC[CH](C)CCC[CH](C)CCC")
    side_chain_unsat = Chem.MolFromSmarts("[CH3]C(C)=CCC[CH](C)=CCC[CH](C)=CCC")
    
    if not (mol.HasSubstructMatch(side_chain_pattern) or mol.HasSubstructMatch(side_chain_unsat)):
        return False, "Missing required isoprenoid chain"
    
    # Check total number of carbons to ensure proper chain length
    carbon_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6])
    if carbon_count < 20:  # Minimum carbons for a tocol
        return False, "Insufficient carbons for tocol structure"
        
    # Additional check for common modifications that don't invalidate tocol status
    valid_modifications = [
        ("acetate", "CC(=O)O"),
        ("succinate", "O=C(O)CCC(=O)O"),
        ("carboxyl", "C(=O)O"),
    ]
    
    # Look for common modifications
    modifications = []
    for mod_name, mod_smarts in valid_modifications:
        pattern = Chem.MolFromSmarts(mod_smarts)
        if mol.HasSubstructMatch(pattern):
            modifications.append(mod_name)
    
    base_message = "Contains chroman-6-ol core with appropriate isoprenoid chain"
    if modifications:
        return True, f"{base_message} (modified with: {', '.join(modifications)})"
    
    return True, base_message


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:39437',
                          'name': 'tocol',
                          'definition': 'A chromanol with a chroman-6-ol '
                                        'skeleton that is substituted at '
                                        'position 2 by a saturated or '
                                        'triply-unsaturated hydrocarbon chain '
                                        'consisting of three isoprenoid units.',
                          'parents': ['CHEBI:23229'],
                          'xrefs': ['CAS:119-98-2', 'Reaxys:1436460'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'CC(=O)O[C@H]1CC[C@]23C[C@@H]1OO[C@@]2(C)C(=O)CCC3(C)C',
                                     'name': 'Talaperoxide B',
                                     'reason': 'No chroman core found'},
                                 {   'smiles': 'C1=CC=C(C(=C1)C#CC2=CC=C(C=C2)[C@H]3[C@H](N[C@@H]3C#N)CO)F',
                                     'name': '(2S,3R,4S)-3-[4-[2-(2-fluorophenyl)ethynyl]phenyl]-4-(hydroxymethyl)-2-azetidinecarbonitrile',
                                     'reason': 'No chroman core found'},
                                 {   'smiles': 'OC(=O)Cc1cn(nc1-c1ccc(Cl)cc1)-c1ccccc1',
                                     'name': 'lonazolac',
                                     'reason': 'No chroman core found'},
                                 {   'smiles': 'C(CCCCCCCC(=O)O)CCCCCCCCC(=O)O',
                                     'name': 'octadecanedioic acid',
                                     'reason': 'No chroman core found'},
                                 {   'smiles': 'O1[C@@H]([C@@H](O)[C@H](O)[C@@H](NC(=O)C)[C@@H]1OC[C@H]2OC(O)[C@H](O)[C@@H](O)[C@H]2O)CO[C@@H]3O[C@@H]([C@@H](O)[C@H](O)[C@H]3NC(=O)C)CO',
                                     'name': 'N-[(2R,3R,4R,5S,6R)-2-[[(2R,3S,4R,5R,6R)-5-Acetamido-3,4-dihydroxy-6-[[(2R,3R,4S,5R)-3,4,5,6-tetrahydroxyoxan-2-yl]methoxy]oxan-2-yl]methoxy]-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]acetamide',
                                     'reason': 'No chroman core found'},
                                 {   'smiles': 'C\\C=C(/C)CC\\C=C(/C)C(=O)OCC(C)(C)CC1=C(O)C(=O)c2ccccc2C1=O',
                                     'name': 'rhinacanthin C',
                                     'reason': 'No chroman core found'},
                                 {   'smiles': 'COC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC(=C(C(=C3)Br)O)Br)OC',
                                     'name': '2,6-dibromo-4-[(6,7-dimethoxy-4-quinazolinyl)amino]phenol',
                                     'reason': 'No chroman core found'},
                                 {   'smiles': 'S(=O)(CC1=CC=CC=C1)C',
                                     'name': 'Methyl benzyl sulfoxide',
                                     'reason': 'No chroman core found'},
                                 {   'smiles': 'C=CCOC1=NS(=O)(=O)c2ccccc12',
                                     'name': 'probenazole',
                                     'reason': 'No chroman core found'},
                                 {   'smiles': 'C(=O)([C@@H](NC(=O)CC)CCSC)[O-]',
                                     'name': 'N-propanoyl-L-methioninate',
                                     'reason': 'No chroman core found'}],
    'sample_false_negatives': [   {   'smiles': 'CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC[C@]1(C)CCc2cc(O)ccc2O1',
                                      'name': 'desmethyl tocotrienol',
                                      'reason': 'Missing required isoprenoid '
                                                'chain'},
                                  {   'smiles': 'C/C(/CC/C=C(\\C)/C(O)=O)=C\\CC/C(/C)=C/CC[C@]1(C)CCC=2C=C(O)C(C)=C(C)C2O1',
                                      'name': "13'-carboxy-gamma-tocotrienol",
                                      'reason': 'Missing required isoprenoid '
                                                'chain'},
                                  {   'smiles': '[H][C@@]1(CC\\C=C(/C)CC\\C=C(/C)CCC=C(C)C)CCc2cc(O)ccc2O1',
                                      'name': 'didesmethyl tocotrienol',
                                      'reason': 'Missing required isoprenoid '
                                                'chain'},
                                  {   'smiles': 'C/C(/CC/C=C(\\C)/C(O)=O)=C/CC/C(/C)=C\\CC[C@]1(C)CCC=2C(C)=C(O)C(C)=C(C)C2O1',
                                      'name': "13'-carboxy-alpha-tocotrienol",
                                      'reason': 'Missing required isoprenoid '
                                                'chain'},
                                  {   'smiles': 'CC(C)=CCC\\C(=C/CC\\C(C)=C\\CC[C@@]1(C)Oc2c(C)cc(O)cc2C=C1)C(O)=O',
                                      'name': '(R)-Sargachromenol',
                                      'reason': 'No chroman core found'},
                                  {   'smiles': 'CC(C)CCCC(C)CCCC(C)CCCC1(C)Oc2c(C)c(C)c(O)cc2C=C1',
                                      'name': 'dehydro-gamma-tocopherol',
                                      'reason': 'No chroman core found'},
                                  {   'smiles': 'CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC[C@]1(C)CCc2c(C)c(O)cc(C)c2O1',
                                      'name': 'beta-tocotrienol',
                                      'reason': 'Missing required isoprenoid '
                                                'chain'},
                                  {   'smiles': 'O1[C@@](CCC2=C1C(=CC(O)=C2)C)(CC/C=C(/CC/C=C(/CC/C=C(\\C)/C(O)=O)\\C)\\C)C',
                                      'name': '(2e,6e,10e)-13-[(2r)-6-hydroxy-2,8-dimethyl-3,4-dihydrochromen-2-yl]-2,6,10-trimethyltrideca-2,6,10-trienoic '
                                              'acid',
                                      'reason': 'Missing required isoprenoid '
                                                'chain'},
                                  {   'smiles': 'CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC[C@]1(C)CCc2c(C)c(O)c(C)c(C)c2O1',
                                      'name': 'alpha-tocotrienol',
                                      'reason': 'Missing required isoprenoid '
                                                'chain'},
                                  {   'smiles': 'O1C(CCC=2C1=C(C(=C(O)C2)C)C)(CCC=C(CCC=C(CCC=C(C)C)C)C)C',
                                      'name': 'I(3)-Tocotrienol',
                                      'reason': 'Missing required isoprenoid '
                                                'chain'}],
    'attempt': 0,
    'reasoning': 'Let me create a program to classify tocols. Based on the '
                 "definition and examples, I'll identify key characteristics:\n"
                 '\n'
                 '1. Base structure is a chromanol (chroman with a hydroxyl '
                 'group)\n'
                 '2. Specifically a chroman-6-ol skeleton\n'
                 '3. Position 2 has a hydrocarbon chain made of 3 isoprenoid '
                 'units\n'
                 '4. Can be saturated or unsaturated\n'
                 '5. May have various methyl substitutions on the aromatic '
                 'ring\n'
                 '6. May have various modifications like acetate or succinate '
                 'esters\n'
                 '\n'
                 "Here's the implementation:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 21,
    'num_false_positives': 4,
    'num_true_negatives': 142262,
    'num_false_negatives': 13,
    'num_negatives': None,
    'precision': 0.84,
    'recall': 0.6176470588235294,
    'f1': 0.711864406779661,
    'accuracy': 0.9998805340829234,
    'negative_predictive_value': 0.9999086276577052}