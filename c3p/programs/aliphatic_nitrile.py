"""
Classifies: CHEBI:80291 aliphatic nitrile
"""
"""
Classifies: CHEBI:22677 aliphatic nitrile
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aliphatic_nitrile(smiles: str):
    """
    Determines if a molecule is an aliphatic nitrile based on its SMILES string.
    An aliphatic nitrile is any nitrile derived from an aliphatic compound.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic nitrile, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for cyano group (C≡N) - more specific pattern that ensures proper connectivity
    cyano_pattern = Chem.MolFromSmarts("[CX2]#[NX1]")
    if cyano_pattern is None:
        return False, "Error in SMARTS pattern"
        
    cyano_matches = mol.GetSubstructMatches(cyano_pattern)
    
    if not cyano_matches:
        return False, "No cyano (C≡N) group found"
    
    # For each cyano group, check if it's attached to an aliphatic carbon
    for match in cyano_matches:
        nitrile_c = mol.GetAtomWithIdx(match[0])  # Get the carbon of C≡N
        
        # Get the carbon atom attached to the nitrile carbon
        for neighbor in nitrile_c.GetNeighbors():
            if neighbor.GetAtomicNum() != 7:  # Skip the nitrogen atom of C≡N
                # Check if this carbon is NOT part of an aromatic system
                if not neighbor.GetIsAromatic():
                    # Check if the parent carbon is not in a ring or is in an aliphatic ring
                    ring_info = mol.GetRingInfo()
                    if not ring_info.IsAtomInRingOfSize(neighbor.GetIdx(), 6) or not neighbor.GetIsAromatic():
                        # Additional check to ensure the whole connected system is not aromatic
                        aromatic_system = False
                        for next_neighbor in neighbor.GetNeighbors():
                            if next_neighbor.GetIsAromatic():
                                aromatic_system = True
                                break
                        
                        if not aromatic_system:
                            return True, "Contains cyano group (C≡N) attached to aliphatic carbon"
    
    return False, "No cyano groups attached to aliphatic carbons found"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:22677',
        'name': 'aliphatic nitrile',
        'definition': 'Any nitrile derived from an aliphatic compound.',
        'parents': ['CHEBI:33577']
    }
}


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:80291',
                          'name': 'aliphatic nitrile',
                          'definition': 'Any nitrile derived from an aliphatic '
                                        'compound.',
                          'parents': ['CHEBI:18379', 'CHEBI:33653'],
                          'xrefs': ['KEGG:C16072'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Error: Python argument types in\n'
               '    Mol.GetSubstructMatches(Mol, NoneType)\n'
               'did not match C++ signature:\n'
               '    GetSubstructMatches(RDKit::ROMol self, RDKit::MolBundle '
               'query, RDKit::SubstructMatchParameters params)\n'
               '    GetSubstructMatches(RDKit::ROMol self, RDKit::ROMol query, '
               'RDKit::SubstructMatchParameters params)\n'
               '    GetSubstructMatches(RDKit::ROMol self, RDKit::MolBundle '
               'query, bool uniquify=True, bool useChirality=False, bool '
               'useQueryQueryMatches=False, unsigned int maxMatches=1000)\n'
               '    GetSubstructMatches(RDKit::ROMol self, RDKit::ROMol query, '
               'bool uniquify=True, bool useChirality=False, bool '
               'useQueryQueryMatches=False, unsigned int maxMatches=1000)\n'
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
    'sample_true_negatives': [   {   'smiles': 'CCS(=O)(=O)NC1=CC2=C(C=C1)OC[C@H](N(C[C@H]([C@H](CN(C2=O)C)OC)C)C(=O)COC)C',
                                     'name': 'N-[(4R,7R,8R)-8-methoxy-5-(2-methoxy-1-oxoethyl)-4,7,10-trimethyl-11-oxo-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]ethanesulfonamide',
                                     'reason': 'No cyano (C≡N) group found'},
                                 {   'smiles': 'C1CCC(C1)(C2=CN(N=N2)CC[C@@H]3CC[C@@H]([C@H](O3)CO)NC(=O)C4=CC=CC=C4F)O',
                                     'name': '2-fluoro-N-[(2S,3S,6S)-6-[2-[4-(1-hydroxycyclopentyl)-1-triazolyl]ethyl]-2-(hydroxymethyl)-3-oxanyl]benzamide',
                                     'reason': 'No cyano (C≡N) group found'},
                                 {   'smiles': 'S(=O)(=O)(NC)CC1=CC2=C(NC=C2CC(O)=O)C=C1',
                                     'name': '1H-Indole-3-acetic acid, '
                                             '5-[[(methylamino)sulfonyl]methyl]-',
                                     'reason': 'No cyano (C≡N) group found'},
                                 {   'smiles': 'O([C@@H]1O[C@@H]([C@H](O)[C@H](O[C@H]2O[C@@H]([C@H](O)[C@H](O)[C@H]2NC(=O)C)CO)[C@H]1O[C@@H]3O[C@H]([C@@H](O)[C@@H](O)[C@@H]3O)C)CO)[C@H]4[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]4CO)OC[C@H]5O[C@H](O)[C@H](NC(=O)C)[C@@H](O[C@@H]6O[C@@H]([C@H](O)[C@H](O)[C@H]6O[C@@H]7O[C@H]([C@@H](O)[C@@H](O)[C@@H]7O)C)CO)[C@H]5O',
                                     'name': 'N-[(2R,3R,4R,5R,6R)-2-[(2S,3R,4S,5S,6R)-2-[(2R,3S,4R,5R,6R)-5-Acetamido-6-[[(2R,3R,4R,5R,6S)-5-acetamido-4-[(2R,3R,4S,5R,6R)-4,5-dihydroxy-6-(hydroxymethyl)-3-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-2-yl]oxy-3,6-dihydroxyoxan-2-yl]methoxy]-4-hydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-5-hydroxy-6-(hydroxymethyl)-3-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-4-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]acetamide',
                                     'reason': 'No cyano (C≡N) group found'},
                                 {   'smiles': 'O([C@H]1[C@H](O)[C@H](O[C@@H](O[C@H]2[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]2CO)O[C@@H]([C@H](O)[C@@H](NC(=O)C)CO)[C@H](O)CO)[C@H]1O)CO[C@H]3O[C@@H]([C@@H](O)[C@H](O[C@H]4O[C@@H]([C@@H](O)[C@H](O)[C@@H]4O[C@H]5O[C@@H]([C@@H](O)[C@H](O)[C@@H]5O)CO)CO)[C@@H]3O)CO[C@H]6O[C@@H]([C@@H](O)[C@H](O)[C@@H]6O[C@H]7O[C@@H]([C@@H](O)[C@H](O)[C@@H]7O)CO)CO)[C@H]8O[C@@H]([C@@H](O)[C@H](O)[C@@H]8O[C@H]9O[C@@H]([C@@H](O)[C@H](O)[C@@H]9O[C@H]%10O[C@@H]([C@@H](O)[C@H](O)[C@@H]%10O)CO)CO)CO',
                                     'name': 'N-[(2S,3R,4R,5S,6R)-2-[(2R,3S,4R,5S)-5-Acetamido-1,2,4,6-tetrahydroxyhexan-3-yl]oxy-5-[(2S,3S,4S,5R,6R)-4-[(2R,3S,4S,5S,6R)-3-[(2S,3S,4S,5S,6R)-4,5-dihydroxy-6-(hydroxymethyl)-3-[(2R,3S,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-6-[[(2S,3S,4S,5R,6R)-4-[(2S,3S,4S,5S,6R)-4,5-dihydroxy-6-(hydroxymethyl)-3-[(2R,3S,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxy-6-[[(2S,3S,4S,5S,6R)-4,5-dihydroxy-6-(hydroxymethyl)-3-[(2R,3S,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxymethyl]-3,5-dihydroxyoxan-2-yl]oxymethyl]-3,5-dihydroxyoxan-2-yl]oxy-4-hydroxy-6-(hydroxymethyl)oxan-3-yl]acetamide',
                                     'reason': 'No cyano (C≡N) group found'},
                                 {   'smiles': 'P(=O)(O[C@H]([C@](O)(C=C[C@@H]1OC(=O)C=C[C@@H]1CC)CCN)C[C@H](O)C=CC=C[C@H]2C[C@@H](OC(=O)C3CCCCC3)CCC2)(O)O',
                                     'name': 'Phoslactomycin E',
                                     'reason': 'No cyano (C≡N) group found'},
                                 {   'smiles': 'ClC1=C(Cl)C=CC=C1C2=N[N+](O)=C(N)N=C2N',
                                     'name': '3,5-diamino-6-(2,3-dichlorophenyl)-2-hydroxy-1,2,4-triazin-2-ium',
                                     'reason': 'No cyano (C≡N) group found'},
                                 {   'smiles': 'P(CC[N+](C)(C)C)(OCC(OC(=O)C)COC(=O)C)(=O)[O-]',
                                     'name': 'diacylglycerol '
                                             '2-trimethylaminoethylphosphonate',
                                     'reason': 'No cyano (C≡N) group found'},
                                 {   'smiles': 'C[C@@H]1CN([C@@H](COC2=C(C=CC(=C2)NC(=O)COC)C(=O)N(C[C@@H]1OC)C)C)CC3=CC(=CC=C3)F',
                                     'name': 'N-[(5R,6R,9R)-8-[(3-fluorophenyl)methyl]-5-methoxy-3,6,9-trimethyl-2-oxo-11-oxa-3,8-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]-2-methoxyacetamide',
                                     'reason': 'No cyano (C≡N) group found'},
                                 {   'smiles': 'NCC(=O)N1[C@H](C(=O)O)C[C@H](C1)O',
                                     'name': 'Gly-Hyp',
                                     'reason': 'No cyano (C≡N) group found'}],
    'sample_false_negatives': [],
    'attempt': 1,
    'reasoning': 'Let me analyze the previous attempt and suggest '
                 'improvements:\n'
                 '\n'
                 '1. The error message indicates that the SMARTS pattern for '
                 "the cyano group wasn't properly created (MolFromSmarts "
                 'returned None).\n'
                 '\n'
                 '2. The previous code only checked if the carbon attached to '
                 "the nitrile was aromatic, but didn't properly verify if it "
                 'was truly aliphatic. An aliphatic carbon should be sp3 or '
                 'sp2 hybridized and not part of any aromatic system.\n'
                 '\n'
                 '3. We should also handle edge cases better and be more '
                 'specific in our pattern matching.\n'
                 '\n'
                 "Here's an improved version:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 31,
    'num_false_positives': 100,
    'num_true_negatives': 16605,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.2366412213740458,
    'recall': 1.0,
    'f1': 0.38271604938271603,
    'accuracy': 0.9940248565965584,
    'negative_predictive_value': 1.0}