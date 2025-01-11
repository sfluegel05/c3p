"""
Classifies: CHEBI:47923 tripeptide
"""
from rdkit import Chem

def is_tripeptide(smiles: str):
    """
    Determines if a molecule is a tripeptide based on its SMILES string.
    A tripeptide consists of three amino acids linked by peptide bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tripeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for peptide bond pattern: N-C(=O)
    peptide_bond_pattern = Chem.MolFromSmarts("N-C(=O)")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    
    # Ensure exactly two peptide bonds for linear tripeptides
    if len(peptide_bond_matches) == 2:
        return True, "Contains exactly 2 peptide bonds typical of a linear tripeptide"

    # Determine if the molecule is cyclic
    kekulized_mol = Chem.Mol(mol)
    contains_ring = kekulized_mol.GetRingInfo().NumRings() > 0

    # If it's cyclic, we may need different logic, but we'll assume cyclic tripeptides for simplicity
    if contains_ring and len(peptide_bond_matches) >= 1:
        return True, f"Contains cyclic structure with peptide bonds, indicative of a cyclic tripeptide"

    # Default to false if none of the expected patterns is properly detected
    return False, f"Found {len(peptide_bond_matches)} peptide bonds, does not match expected tripeptide structure"

# Example usage
# result, reason = is_tripeptide("CC[C@H](C)[C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)N[C@@H](CO)C(O)=O")
# print(result, reason)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:47923',
                          'name': 'tripeptide',
                          'definition': 'Any oligopeptide that consists of '
                                        'three amino-acid residues connected '
                                        'by peptide linkages.',
                          'parents': ['CHEBI:25676'],
                          'xrefs': ['KEGG:C00316', 'Wikipedia:Tripeptide'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               "Error: '>' not supported between instances of "
               "'_vectNSt3__16vectorIiNS_9allocatorIiEEEE' and 'int'\n"
               'Attempt failed: F1 score of 0 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: NONE\n'
               'False positives: NONE\n'
               'False negatives: NONE\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'P(OCC(OC(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCC/C=C\\C/C=C\\C/C=C\\CCCCCCCC)(OCCN(C)C)(O)=O',
                                     'name': 'Pe-nme2(20:3(5Z,8Z,11Z)/22:6(4Z,7Z,10Z,13Z,16Z,19Z))',
                                     'reason': 'Found 0 peptide bonds, does '
                                               'not match expected tripeptide '
                                               'structure'},
                                 {   'smiles': 'O1[C@]2([C@@H](O)[C@H](O)[C@@H](O)[C@@]1(OCC[C@H](CCC=C(C(OC[C@]3(O)[C@@H](O)[C@@](OC3)(OC2)[H])=O)C)C)[H])[H]',
                                     'name': 'Urceolide',
                                     'reason': 'Found 0 peptide bonds, does '
                                               'not match expected tripeptide '
                                               'structure'},
                                 {   'smiles': 'CC1=CC(=C(C=C1)C)C(=O)CSC2=NN=C(S2)C',
                                     'name': '1-(2,5-dimethylphenyl)-2-[(5-methyl-1,3,4-thiadiazol-2-yl)thio]ethanone',
                                     'reason': 'Found 0 peptide bonds, does '
                                               'not match expected tripeptide '
                                               'structure'},
                                 {   'smiles': 'O=C1C=C([C@@]2(C(C[C@@H](C2)O)(C)C)C)CC[C@@]1(O)C',
                                     'name': 'Enokipodin H',
                                     'reason': 'Found 0 peptide bonds, does '
                                               'not match expected tripeptide '
                                               'structure'},
                                 {   'smiles': 'CCCC[C@](Cn1cncn1)(C#N)c1ccc(Cl)cc1',
                                     'name': '(S)-myclobutanil',
                                     'reason': 'Found 0 peptide bonds, does '
                                               'not match expected tripeptide '
                                               'structure'},
                                 {   'smiles': 'OC(=O)CCCCCCC#CCCCC',
                                     'name': '8-tridecynoic acid',
                                     'reason': 'Found 0 peptide bonds, does '
                                               'not match expected tripeptide '
                                               'structure'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@H](O[C@@H]3[C@@H](CO)O[C@H](O[C@@H]4[C@@H](CO)O[C@H](O)[C@H](O)[C@H]4O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O',
                                     'name': 'alpha-maltotetraose',
                                     'reason': 'Found 0 peptide bonds, does '
                                               'not match expected tripeptide '
                                               'structure'},
                                 {   'smiles': 'SC[C@H](N)C(=O)N[C@H](C(=O)N[C@@H](CCC(=O)N)C(O)=O)CO',
                                     'name': 'Cys-Ser-Gln',
                                     'reason': 'Found 3 peptide bonds, does '
                                               'not match expected tripeptide '
                                               'structure'},
                                 {   'smiles': 'S(C(C)C(O)=O)CC1=CC=CC=C1',
                                     'name': '2-(Benzylthio)propanoic acid',
                                     'reason': 'Found 0 peptide bonds, does '
                                               'not match expected tripeptide '
                                               'structure'},
                                 {   'smiles': '[H][C@@]1(CNc2nc(N)[nH]c(=O)c2N1)[C@@H](O)[C@H](C)O',
                                     'name': 'sapropterin',
                                     'reason': 'Found 0 peptide bonds, does '
                                               'not match expected tripeptide '
                                               'structure'}],
    'sample_false_negatives': [   {   'smiles': 'CC(C)C[C@H](N)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](CC(N)=O)C(O)=O',
                                      'name': 'Leu-Leu-Asn',
                                      'reason': 'Found 3 peptide bonds, does '
                                                'not match expected tripeptide '
                                                'structure'},
                                  {   'smiles': 'C(=O)([C@@H](N)CCSC)N[C@H](C(=O)N[C@H](C(=O)O)CO)CC(=O)N',
                                      'name': 'Met-Asn-Ser',
                                      'reason': 'Found 3 peptide bonds, does '
                                                'not match expected tripeptide '
                                                'structure'},
                                  {   'smiles': 'NC(=O)CC[C@@H](C(=O)N[C@H](C(N[C@H](C(=O)O)CCCNC(=N)N)=O)CCCCN)N',
                                      'name': 'Gln-Lys-Arg',
                                      'reason': 'Found 3 peptide bonds, does '
                                                'not match expected tripeptide '
                                                'structure'},
                                  {   'smiles': 'CC(=O)N[C@@H](CC(O)=O)C(=O)N[C@@H](CCC(O)=O)C(=O)N[C@@H](CCC(O)=O)C(O)=O',
                                      'name': 'Ac-Asp-Glu-Glu',
                                      'reason': 'Found 3 peptide bonds, does '
                                                'not match expected tripeptide '
                                                'structure'},
                                  {   'smiles': 'CCCC[C@H](NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC(C)C)NC(C)=O)C=O',
                                      'name': 'acetylleucyl-leucyl-norleucinal',
                                      'reason': 'Found 3 peptide bonds, does '
                                                'not match expected tripeptide '
                                                'structure'},
                                  {   'smiles': 'C(=O)([C@@H](N)CC(C)C)N[C@H](C(=O)N[C@H](C(=O)O)CO)CC(=O)N',
                                      'name': 'Leu-Asn-Ser',
                                      'reason': 'Found 3 peptide bonds, does '
                                                'not match expected tripeptide '
                                                'structure'},
                                  {   'smiles': 'CC(C)C[C@H](N)C(=O)N[C@@H]([C@@H](C)O)C(=O)N[C@@H](CCC(N)=O)C(O)=O',
                                      'name': 'Leu-Thr-Gln',
                                      'reason': 'Found 3 peptide bonds, does '
                                                'not match expected tripeptide '
                                                'structure'},
                                  {   'smiles': 'C([C@@H](NC([C@@H](NC([C@@H](NC(C)=O)CC(C)C)=O)CC(C)C)=O)CCCNC(=N)N)(=O)[H]',
                                      'name': 'leupeptin',
                                      'reason': 'Found 3 peptide bonds, does '
                                                'not match expected tripeptide '
                                                'structure'},
                                  {   'smiles': 'N[C@@H](CCC(=O)N[C@@H](CS)C(=O)NCC(N)=O)C(O)=O',
                                      'name': 'glutathione amide',
                                      'reason': 'Found 3 peptide bonds, does '
                                                'not match expected tripeptide '
                                                'structure'},
                                  {   'smiles': 'C(=O)([C@@H](N)CCC(=O)O)N[C@H](C(=O)N[C@H](C(=O)O)CCC(N)=O)CCC(=O)O',
                                      'name': 'Glu-Glu-Gln',
                                      'reason': 'Found 3 peptide bonds, does '
                                                'not match expected tripeptide '
                                                'structure'}],
    'attempt': 2,
    'reasoning': 'The previous program attempted to classify tripeptides by '
                 'identifying peptide bond patterns and counting these bonds '
                 'to determine if a molecule fits the definition of a '
                 'tripeptide, which generally consists of three amino acids '
                 'linked by peptide bonds. The error appears to stem from '
                 'using a method `GetSSSR`, resulting in an unsupported '
                 'comparison further downstream. Additionally, there was a '
                 'lack of sufficient checks for specific tripeptide patterns, '
                 'like determining the presence of amino acids count and '
                 'structure patterns.\n'
                 '\n'
                 "To improve upon this, we'll focus on:\n"
                 '1. Correcting the use of the ring detection method.\n'
                 "2. Ensuring we're detecting exactly three amino acids "
                 'connected by peptide bonds.\n'
                 '3. Improving pattern matching for both linear and cyclic '
                 'tripeptides.\n'
                 '\n'
                 "Let's now redefine the function to address these issues:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 155,
    'num_false_positives': 100,
    'num_true_negatives': 150,
    'num_false_negatives': 11,
    'num_negatives': None,
    'precision': 0.6078431372549019,
    'recall': 0.9337349397590361,
    'f1': 0.7363420427553443,
    'accuracy': 0.7331730769230769,
    'negative_predictive_value': 0.9316770186335404}