"""
Classifies: CHEBI:29348 fatty amide
"""
"""
Classifies: fatty amide compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_fatty_amide(smiles: str):
    """
    Determines if a molecule is a fatty amide based on its SMILES string.
    A fatty amide is a monocarboxylic acid amide derived from a fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty amide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for amide group
    amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3]")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    
    if not amide_matches:
        return False, "No amide group found"
    
    # Get number of amide groups
    num_amides = len(amide_matches)
    
    # Check carbon chain length
    carbon_chain = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    chain_matches = mol.GetSubstructMatches(carbon_chain)
    
    if not chain_matches:
        return False, "No significant carbon chain found"
    
    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 4:
        return False, "Carbon chain too short for fatty amide"
        
    # Count rotatable bonds to check chain flexibility
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    
    # Additional checks for specific cases that should be excluded
    
    # Check for cyclic peptides (multiple amides in cycle)
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        amide_in_ring = 0
        for amide_match in amide_matches:
            if amide_match[0] in ring:
                amide_in_ring += 1
        if amide_in_ring > 1:
            return False, "Appears to be a cyclic peptide"
    
    # Check for peptides (multiple amides in sequence)
    peptide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[CX4][NX3][CX3](=[OX1])")
    if mol.HasSubstructMatch(peptide_pattern):
        return False, "Appears to be a peptide"
        
    # Positive classification criteria
    if num_amides == 1:
        if c_count >= 4 and n_rotatable >= 2:
            return True, "Contains single amide group with appropriate carbon chain"
    else:
        # For molecules with multiple amides, check if they're separated by long chains
        if c_count/num_amides >= 4 and n_rotatable >= num_amides*2:
            return True, "Contains well-separated amide groups with appropriate carbon chains"
            
    return False, "Does not match fatty amide criteria"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:29348',
                          'name': 'fatty amide',
                          'definition': 'A monocarboxylic acid amide derived '
                                        'from a fatty acid.',
                          'parents': ['CHEBI:29347', 'CHEBI:61697'],
                          'xrefs': ['KEGG:C02244', 'LIPID_MAPS_class:LMFA08'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'Cc1cc(C)c2Cc3c(C)cc(C)c(Cc4c(C)cc(C)c(Cc5c(C)cc(C)c(Cc1c2C)c5C)c4C)c3C',
                                     'name': '4,6,10,12,16,18,22,24,25,26,27,28-dodecamethylcalix[4]arene',
                                     'reason': 'No amide group found'},
                                 {   'smiles': 'O=C1C(=C2C=C3[C@]([C@@H](C(C)C)[C@@H]([C@H]3O)OC(=O)C)(C)CC[C@]2(C)CC1)COC(=O)C',
                                     'name': 'Dahliane E',
                                     'reason': 'No amide group found'},
                                 {   'smiles': 'O[C@@H]([C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)N[C@@H](CC(C)C)C(O)=O)C',
                                     'name': 'Glu-Thr-Leu',
                                     'reason': 'Appears to be a peptide'},
                                 {   'smiles': 'CCOc1ccc(NC(=O)C(C)O)cc1',
                                     'name': 'p-Lactophenetide',
                                     'reason': 'No significant carbon chain '
                                               'found'},
                                 {   'smiles': 'O=C(OC1C(O)C(OC(C1O)CO)OCC2OC(OCC(OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)C(O)C(C2O)O)CCCCCCCCCCCCCCC',
                                     'name': '[(2S)-2-hexadecanoyloxy-3-[(2S,3R,4S,5S,6R)-6-[[(2S,3R,4S,5S,6R)-4-hexadecanoyloxy-3,5-dihydroxy-6-(hydroxymethyl)tetrahydropyran-2-yl]oxymethyl]-3,4,5-trihydroxy-tetrahydropyran-2-yl]oxy-propyl] '
                                             '(9E,12E,15E)-octadeca-9,12,15-trienoate',
                                     'reason': 'No amide group found'},
                                 {   'smiles': 'O=C1C2=C(O)C=C(OC)C=C2C(=O)C3=C1[C@@H]([C@@H](O)[C@]([C@@H]3O)(O)C)C',
                                     'name': 'Altersolanol G',
                                     'reason': 'No amide group found'},
                                 {   'smiles': 'OC(C(O)C/C=C\\C/C=C\\C/C=C\\CC)C/C=C\\C/C=C\\CCC(O)=O',
                                     'name': '10,11-DiHDPE',
                                     'reason': 'No amide group found'},
                                 {   'smiles': '[Na+].[H][C@]12SCC(C)=C(N1C(=O)[C@H]2NC(=O)[C@H](N)c1ccccc1)C([O-])=O',
                                     'name': 'cephalexin sodium',
                                     'reason': 'Appears to be a peptide'},
                                 {   'smiles': 'C1=CC=C(C=C1)COC2=CC=CC=C2C=NN3C(=NN=N3)N',
                                     'name': '1-[(2-phenylmethoxyphenyl)methylideneamino]-5-tetrazolamine',
                                     'reason': 'No amide group found'},
                                 {   'smiles': 'O[C@]1(C([C@H]2[C@](C3=C([C@]4([C@]([C@@H]([C@@H]([C@@H](O)CC=C(C)C)C)CC4)(C)CC3)C)CC2)(C)CC1)(C)C)C',
                                     'name': 'Astraeusin N',
                                     'reason': 'No amide group found'}],
    'sample_false_negatives': [   {   'smiles': 'CCC(N)=O',
                                      'name': 'propionamide',
                                      'reason': 'No significant carbon chain '
                                                'found'},
                                  {   'smiles': 'CC(=O)NCCC1=CC(=C(C=C1)O)O',
                                      'name': 'N-acetyldopamine',
                                      'reason': 'No significant carbon chain '
                                                'found'},
                                  {   'smiles': 'C1=CC(=CC=2C(=CNC12)CCNC(C)=O)O',
                                      'name': 'N-acetylserotonin',
                                      'reason': 'No significant carbon chain '
                                                'found'},
                                  {   'smiles': 'O=C(N)CCCCC(=O)N',
                                      'name': 'Adipamide',
                                      'reason': 'Does not match fatty amide '
                                                'criteria'},
                                  {   'smiles': 'O(C(C)(C)C)C(=O)NC(C)C(=O)C(=O)N',
                                      'name': 'Sakacin P',
                                      'reason': 'Does not match fatty amide '
                                                'criteria'},
                                  {   'smiles': 'O=C(NCCC1=CC=CC=C1)[C@@H](N(C(=O)C(=O)C(CC)C)C)C(C)C',
                                      'name': 'Nevaltophin A',
                                      'reason': 'Appears to be a peptide'},
                                  {   'smiles': 'O=C(N([C@H](C(=O)NCCC1=CC=CC=C1)C(C)C)C)C(=O)C(C)C',
                                      'name': 'Nevaltophin B',
                                      'reason': 'Appears to be a peptide'},
                                  {   'smiles': 'CC(CCCCCCCCCCCC([C@@H]([NH3+])CS(=O)(=O)[O-])=O)C',
                                      'name': '2-amino-3-oxo-15-methylhexadecane-1-sulfonate',
                                      'reason': 'No amide group found'},
                                  {   'smiles': 'CCCCCCCCCCCCCCCC([C@@H]([NH3+])CS(=O)(=O)[O-])=O',
                                      'name': '2-amino-3-oxo-octadecane-1-sulfonate',
                                      'reason': 'No amide group found'},
                                  {   'smiles': 'CC(CCCCCCCCCCCC([C@@H]([NH3+])CS(=O)(=O)[O-])O)C',
                                      'name': '2-amino-3-hydroxy-15-methylhexadecane-1-sulfonate',
                                      'reason': 'No amide group found'}],
    'attempt': 0,
    'reasoning': "I'll create a program to classify fatty amides. Let's break "
                 'this down:\n'
                 '\n'
                 '1. A fatty amide must have:\n'
                 '- An amide group (-C(=O)-N)\n'
                 '- A fatty acid-derived chain (typically >4 carbons)\n'
                 '\n'
                 '2. Key patterns to check:\n'
                 '- Amide group presence\n'
                 '- Chain length\n'
                 '- Absence of conflicting functional groups\n'
                 '- Presence of typical fatty acid features (long alkyl '
                 'chains, possible unsaturation)\n'
                 '\n'
                 "Here's the implementation:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 294,
    'num_false_positives': 100,
    'num_true_negatives': 305,
    'num_false_negatives': 10,
    'num_negatives': None,
    'precision': 0.7461928934010152,
    'recall': 0.9671052631578947,
    'f1': 0.8424068767908309,
    'accuracy': 0.844851904090268,
    'negative_predictive_value': 0.9682539682539683}