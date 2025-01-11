"""
Classifies: CHEBI:32955 epoxide
"""
"""
Classifies: CHEBI:35762 epoxide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_epoxide(smiles: str):
    """
    Determines if a molecule contains an epoxide group based on its SMILES string.
    An epoxide is a cyclic ether where the oxygen atom forms part of a 3-membered ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains an epoxide group, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Try to sanitize the molecule
    try:
        Chem.SanitizeMol(mol)
    except:
        return False, "Could not sanitize molecule"

    # SMARTS pattern for epoxide - a three-membered ring containing an oxygen
    # [OX2r3] - oxygen with 2 connections in a 3-membered ring
    # 1 - first connection point
    # [#6r3] - any carbon in a 3-membered ring
    # [#6r3] - second carbon in a 3-membered ring
    # The whole pattern must form a 3-membered ring
    epoxide_pattern = Chem.MolFromSmarts('[OX2r3]1[#6r3][#6r3]1')
    
    if epoxide_pattern is None:
        return False, "Error in SMARTS pattern"
    
    # Look for matches with stereochemistry consideration
    matches = mol.GetSubstructMatches(epoxide_pattern, useChirality=True)
    
    if not matches:
        return False, "No epoxide group found"
        
    # Verify each match is actually a 3-membered ring containing an oxygen
    ring_info = mol.GetRingInfo()
    for match in matches:
        # Get the ring size for these atoms
        rings = ring_info.AtomRings()
        for ring in rings:
            if all(idx in ring for idx in match):
                if len(ring) == 3:  # Confirm it's a 3-membered ring
                    # Verify one atom is oxygen and others are carbons
                    atoms = [mol.GetAtomWithIdx(idx) for idx in match]
                    atom_types = [atom.GetAtomicNum() for atom in atoms]
                    if atom_types.count(8) == 1:  # One oxygen
                        return True, "Contains epoxide group (3-membered ring with oxygen)"
    
    return False, "No valid epoxide group found"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:35762',
        'name': 'epoxide',
        'definition': 'Any cyclic ether in which the oxygen atom forms part of a 3-membered ring.',
        'parents': ['CHEBI:33641']
    }
}


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:32955',
                          'name': 'epoxide',
                          'definition': 'Any cyclic ether in which the oxygen '
                                        'atom forms part of a 3-membered ring.',
                          'parents': ['CHEBI:37407'],
                          'xrefs': [   'KEGG:C00722',
                                       'PMID:10891060',
                                       'Wikipedia:Epoxide'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Error: Python argument types in\n'
               '    Mol.HasSubstructMatch(Mol, NoneType)\n'
               'did not match C++ signature:\n'
               '    HasSubstructMatch(RDKit::ROMol self, RDKit::MolBundle '
               'query, RDKit::SubstructMatchParameters params=True)\n'
               '    HasSubstructMatch(RDKit::ROMol self, RDKit::ROMol query, '
               'RDKit::SubstructMatchParameters params)\n'
               '    HasSubstructMatch(RDKit::ROMol self, RDKit::MolBundle '
               'query, bool recursionPossible=True, bool useChirality=False, '
               'bool useQueryQueryMatches=False)\n'
               '    HasSubstructMatch(RDKit::ROMol self, RDKit::ROMol query, '
               'bool recursionPossible=True, bool useChirality=False, bool '
               'useQueryQueryMatches=False)\n'
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
    'sample_true_negatives': [   {   'smiles': 'Cc1cc(C)c2Cc3c(C)cc(C)c(Cc4c(C)cc(C)c(Cc5c(C)cc(C)c(Cc1c2C)c5C)c4C)c3C',
                                     'name': '4,6,10,12,16,18,22,24,25,26,27,28-dodecamethylcalix[4]arene',
                                     'reason': 'No epoxide group found'},
                                 {   'smiles': 'O=C1C(=C2C=C3[C@]([C@@H](C(C)C)[C@@H]([C@H]3O)OC(=O)C)(C)CC[C@]2(C)CC1)COC(=O)C',
                                     'name': 'Dahliane E',
                                     'reason': 'No epoxide group found'},
                                 {   'smiles': 'O[C@@H]([C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)N[C@@H](CC(C)C)C(O)=O)C',
                                     'name': 'Glu-Thr-Leu',
                                     'reason': 'No epoxide group found'},
                                 {   'smiles': 'CCOc1ccc(NC(=O)C(C)O)cc1',
                                     'name': 'p-Lactophenetide',
                                     'reason': 'No epoxide group found'},
                                 {   'smiles': 'O=C1NCC=2C1=C3C(N([C@@H]4O[C@H]([C@@H](O)[C@H]([C@H]4OC)O)C)C5=C3C=CC=C5)=C6NC7=C(C26)C=CC=C7',
                                     'name': "3'-epi-5'-methoxy-K252d",
                                     'reason': 'No epoxide group found'},
                                 {   'smiles': 'O=C(OC1C(O)C(OC(C1O)CO)OCC2OC(OCC(OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)C(O)C(C2O)O)CCCCCCCCCCCCCCC',
                                     'name': '[(2S)-2-hexadecanoyloxy-3-[(2S,3R,4S,5S,6R)-6-[[(2S,3R,4S,5S,6R)-4-hexadecanoyloxy-3,5-dihydroxy-6-(hydroxymethyl)tetrahydropyran-2-yl]oxymethyl]-3,4,5-trihydroxy-tetrahydropyran-2-yl]oxy-propyl] '
                                             '(9E,12E,15E)-octadeca-9,12,15-trienoate',
                                     'reason': 'No epoxide group found'},
                                 {   'smiles': 'O=C1C2=C(O)C=C(OC)C=C2C(=O)C3=C1[C@@H]([C@@H](O)[C@]([C@@H]3O)(O)C)C',
                                     'name': 'Altersolanol G',
                                     'reason': 'No epoxide group found'},
                                 {   'smiles': '[H][C@]1(O[C@](O)(C[C@H](O)[C@H]1NC(=O)CO)C([O-])=O)[C@H](O)[C@H](O)CO',
                                     'name': 'N-glycoloyl-alpha-neuraminate',
                                     'reason': 'No epoxide group found'},
                                 {   'smiles': 'OC(C(O)C/C=C\\C/C=C\\C/C=C\\CC)C/C=C\\C/C=C\\CCC(O)=O',
                                     'name': '10,11-DiHDPE',
                                     'reason': 'No epoxide group found'},
                                 {   'smiles': '[Na+].[H][C@]12SCC(C)=C(N1C(=O)[C@H]2NC(=O)[C@H](N)c1ccccc1)C([O-])=O',
                                     'name': 'cephalexin sodium',
                                     'reason': 'No epoxide group found'}],
    'sample_false_negatives': [   {   'smiles': 'C([O-])(=O)CCCCCCC/C=C\\[C@@H]1[C@@H]2[C@@H](C[C@]1([C@H](CC)O2)[H])O',
                                      'name': 'plasmodiophorol B(1-)',
                                      'reason': 'No epoxide group found'},
                                  {   'smiles': 'C1[C@]2([C@]3([C@@]4([C@@]([C@H]([C@@]1(C)C([C@H]4O)=O)S(O)(=O)=O)(O[C@@]2([C@H](O)C3)[H])[H])CO)C)O',
                                      'name': 'DON-sulfonate 1',
                                      'reason': 'No epoxide group found'},
                                  {   'smiles': 'CC(\\C=C\\C=C(/C)C1OC2(C)CC(O)CC(C)(C)C2=C1)=C/C=C/C=C(C)/C=C/C=C(\\C)C1OC2(C)CC(O)CC(C)(C)C2=C1',
                                      'name': 'auroxanthin',
                                      'reason': 'No epoxide group found'}],
    'attempt': 2,
    'reasoning': 'Let me analyze the previous attempt and suggest '
                 'improvements:\n'
                 '\n'
                 '1. Error Analysis:\n'
                 '- The error occurred because one of the SMARTS patterns '
                 'failed to parse correctly (returned None)\n'
                 '- We need to verify each SMARTS pattern is valid before '
                 'using it\n'
                 '- The patterns may not have been comprehensive enough to '
                 'catch all epoxide variations\n'
                 '\n'
                 '2. Issues with Previous Approach:\n'
                 '- The ring detection logic was overly complex and '
                 'potentially redundant\n'
                 '- SMARTS patterns could be simplified and made more robust\n'
                 '- Some stereochemistry considerations might have been '
                 'missing\n'
                 "- The code didn't properly handle error cases\n"
                 '\n'
                 '3. Improvements Needed:\n'
                 '- Simplify the SMARTS pattern to a single, robust pattern\n'
                 '- Add proper error handling for SMARTS parsing\n'
                 '- Consider stereochemistry in matching\n'
                 '- Simplify the overall logic\n'
                 '- Add more validation checks\n'
                 '\n'
                 "Here's the improved code:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 650,
    'num_false_positives': 100,
    'num_true_negatives': 6351,
    'num_false_negatives': 3,
    'num_negatives': None,
    'precision': 0.8666666666666667,
    'recall': 0.9954058192955589,
    'f1': 0.9265858873841769,
    'accuracy': 0.9855011261261262,
    'negative_predictive_value': 0.9995278564683664}