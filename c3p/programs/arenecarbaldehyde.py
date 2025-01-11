"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
"""
Classifies: arenecarbaldehyde
Definition: Any aldehyde in which the carbonyl group is attached to an aromatic moiety.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde based on its SMILES string.
    An arenecarbaldehyde has an aldehyde group (-CHO) directly attached to an aromatic ring.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (is_arenecarbaldehyde, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all aldehyde groups
    # Match both explicit and implicit hydrogen representations
    aldehyde_pattern = Chem.MolFromSmarts("[$([CH](=O))]")
    matches = mol.GetSubstructMatches(aldehyde_pattern)
    
    if not matches:
        return False, "No aldehyde group found"

    # For each aldehyde carbon, check if it's connected to an aromatic system
    for match in matches:
        aldehyde_carbon = mol.GetAtomWithIdx(match[0])
        
        # Get all neighboring atoms of the aldehyde carbon
        for neighbor in aldehyde_carbon.GetNeighbors():
            # Skip the oxygen of the aldehyde group
            if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() == 0:
                continue
                
            # Check if the neighbor is part of an aromatic system
            if neighbor.GetIsAromatic():
                return True, "Contains aldehyde group directly attached to aromatic ring"
            
            # Check if neighbor is part of a ring that might be aromatic
            ring_info = mol.GetRingInfo()
            atom_rings = ring_info.AtomRings()
            
            for ring in atom_rings:
                if neighbor.GetIdx() in ring:
                    # Check if any atom in this ring is aromatic
                    ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
                    if any(atom.GetIsAromatic() for atom in ring_atoms):
                        return True, "Contains aldehyde group directly attached to aromatic ring system"

    return False, "No aldehyde group attached to aromatic ring"

def test_smiles():
    """Test function with known examples"""
    test_cases = [
        ("O=CC1=CC=C(N)C=C1", True),  # p-aminobenzaldehyde
        ("Cc1ccccc1C=O", True),        # o-tolualdehyde
        ("O=CC=1N(CC)C=CC1", True),    # 1-Ethyl-1H-pyrrole-2-carboxaldehyde
        ("CCN(CC)c1ccc(C=O)cc1", True),# 4-(diethylamino)benzaldehyde
        ("CC(=O)CC(=O)CC", False),     # Not an arenecarbaldehyde
        ("c1ccccc1", False),           # Not an aldehyde
    ]
    
    for smiles, expected in test_cases:
        result, reason = is_arenecarbaldehyde(smiles)
        print(f"SMILES: {smiles}")
        print(f"Expected: {expected}, Got: {result}")
        print(f"Reason: {reason}\n")

if __name__ == "__main__":
    test_smiles()


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33855',
                          'name': 'arenecarbaldehyde',
                          'definition': 'Any aldehyde in which the carbonyl '
                                        'group is attached to an aromatic '
                                        'moiety.',
                          'parents': ['CHEBI:17478', 'CHEBI:33659'],
                          'xrefs': ['KEGG:C00193', 'MetaCyc:Aryl-Aldehyde'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               "Error: 'RingInfo' object has no attribute 'IsAtomInRing'\n"
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
                                     'reason': 'No aldehyde group found'},
                                 {   'smiles': 'O=C1C(=C2C=C3[C@]([C@@H](C(C)C)[C@@H]([C@H]3O)OC(=O)C)(C)CC[C@]2(C)CC1)COC(=O)C',
                                     'name': 'Dahliane E',
                                     'reason': 'No aldehyde group found'},
                                 {   'smiles': 'O[C@@H]([C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)N[C@@H](CC(C)C)C(O)=O)C',
                                     'name': 'Glu-Thr-Leu',
                                     'reason': 'No aldehyde group found'},
                                 {   'smiles': 'CCOc1ccc(NC(=O)C(C)O)cc1',
                                     'name': 'p-Lactophenetide',
                                     'reason': 'No aldehyde group found'},
                                 {   'smiles': 'O=C1NCC=2C1=C3C(N([C@@H]4O[C@H]([C@@H](O)[C@H]([C@H]4OC)O)C)C5=C3C=CC=C5)=C6NC7=C(C26)C=CC=C7',
                                     'name': "3'-epi-5'-methoxy-K252d",
                                     'reason': 'No aldehyde group found'},
                                 {   'smiles': 'O=C(OC1C(O)C(OC(C1O)CO)OCC2OC(OCC(OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)C(O)C(C2O)O)CCCCCCCCCCCCCCC',
                                     'name': '[(2S)-2-hexadecanoyloxy-3-[(2S,3R,4S,5S,6R)-6-[[(2S,3R,4S,5S,6R)-4-hexadecanoyloxy-3,5-dihydroxy-6-(hydroxymethyl)tetrahydropyran-2-yl]oxymethyl]-3,4,5-trihydroxy-tetrahydropyran-2-yl]oxy-propyl] '
                                             '(9E,12E,15E)-octadeca-9,12,15-trienoate',
                                     'reason': 'No aldehyde group found'},
                                 {   'smiles': 'O=C1C2=C(O)C=C(OC)C=C2C(=O)C3=C1[C@@H]([C@@H](O)[C@]([C@@H]3O)(O)C)C',
                                     'name': 'Altersolanol G',
                                     'reason': 'No aldehyde group found'},
                                 {   'smiles': '[H][C@]1(O[C@](O)(C[C@H](O)[C@H]1NC(=O)CO)C([O-])=O)[C@H](O)[C@H](O)CO',
                                     'name': 'N-glycoloyl-alpha-neuraminate',
                                     'reason': 'No aldehyde group found'},
                                 {   'smiles': 'OC(C(O)C/C=C\\C/C=C\\C/C=C\\CC)C/C=C\\C/C=C\\CCC(O)=O',
                                     'name': '10,11-DiHDPE',
                                     'reason': 'No aldehyde group found'},
                                 {   'smiles': '[Na+].[H][C@]12SCC(C)=C(N1C(=O)[C@H]2NC(=O)[C@H](N)c1ccccc1)C([O-])=O',
                                     'name': 'cephalexin sodium',
                                     'reason': 'No aldehyde group found'}],
    'sample_false_negatives': [],
    'attempt': 2,
    'reasoning': 'Let me analyze the previous code and propose improvements:\n'
                 '\n'
                 '1. Error Analysis:\n'
                 "- The main error 'RingInfo' object has no attribute "
                 "'IsAtomInRing' suggests we used an incorrect method to check "
                 'if an atom is in a ring\n'
                 '- The correct method is `IsAtomInRingOfSize()` or checking '
                 'the atom rings directly\n'
                 '\n'
                 '2. Potential Issues:\n'
                 '- The aromatic detection might be too restrictive\n'
                 '- We need to handle both explicit and implicit hydrogens '
                 'better\n'
                 '- The ring detection logic could be simplified\n'
                 '\n'
                 '3. Improvements:\n'
                 '- Use simpler SMARTS patterns to detect aldehyde groups\n'
                 '- Simplify the ring and aromaticity checking logic\n'
                 '- Add better handling of different aromatic systems '
                 '(benzene, pyrrole, furan, etc.)\n'
                 "- Use RDKit's built-in aromaticity perception more "
                 'effectively\n'
                 '\n'
                 "Here's the improved code:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 136,
    'num_false_positives': 100,
    'num_true_negatives': 25509,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.576271186440678,
    'recall': 1.0,
    'f1': 0.7311827956989247,
    'accuracy': 0.9961157506311905,
    'negative_predictive_value': 1.0}