"""
Classifies: CHEBI:26888 tetrachlorobenzene
"""
"""
Classifies: tetrachlorobenzene
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tetrachlorobenzene(smiles: str):
    """
    Determines if a molecule contains a tetrachlorobenzene substructure.
    Tetrachlorobenzene is a benzene ring with exactly 4 chlorine atoms attached.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains tetrachlorobenzene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all benzene rings
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    benzene_matches = mol.GetSubstructMatches(benzene_pattern)
    
    if not benzene_matches:
        return False, "No benzene ring found"
    
    # For each benzene ring, check if it has exactly 4 chlorines attached
    for benzene_match in benzene_matches:
        benzene_atoms = set(benzene_match)
        chlorine_count = 0
        
        # Check each carbon in the benzene ring for attached chlorines
        for atom_idx in benzene_match:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 17:  # Chlorine atomic number
                    chlorine_count += 1
        
        if chlorine_count == 4:
            return True, "Found benzene ring with exactly 4 chlorine atoms"
            
    return False, "No benzene ring with exactly 4 chlorines found"

def test_examples():
    """Test function with known examples"""
    examples = [
        "Clc1cc(Cl)c(Cl)cc1Cl",  # 1,2,4,5-tetrachlorobenzene
        "Clc1ccc(Cl)c(Cl)c1Cl",  # 1,2,3,4-tetrachlorobenzene
        "Clc1cc(Cl)c(Cl)c(Cl)c1",  # 1,2,3,5-tetrachlorobenzene
        "CC1=CC=CC=C1",  # Toluene (negative example)
        "Clc1c(Cl)c(Cl)c(Cl)c(Cl)c1Cl"  # Hexachlorobenzene (negative example)
    ]
    
    for smiles in examples:
        result, reason = is_tetrachlorobenzene(smiles)
        print(f"SMILES: {smiles}")
        print(f"Result: {result}")
        print(f"Reason: {reason}\n")

if __name__ == "__main__":
    test_examples()


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26888',
                          'name': 'tetrachlorobenzene',
                          'definition': 'Any member of the class of '
                                        'chlorobenzenes carrying four chloro '
                                        'groups at unspecified positions.',
                          'parents': ['CHEBI:23132'],
                          'xrefs': ['CAS:12408-10-5'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'Cc1cc(C)c2Cc3c(C)cc(C)c(Cc4c(C)cc(C)c(Cc5c(C)cc(C)c(Cc1c2C)c5C)c4C)c3C',
                                     'name': '4,6,10,12,16,18,22,24,25,26,27,28-dodecamethylcalix[4]arene',
                                     'reason': 'No benzene ring with exactly 4 '
                                               'chlorines found'},
                                 {   'smiles': 'O=C1C(=C2C=C3[C@]([C@@H](C(C)C)[C@@H]([C@H]3O)OC(=O)C)(C)CC[C@]2(C)CC1)COC(=O)C',
                                     'name': 'Dahliane E',
                                     'reason': 'No benzene ring found'},
                                 {   'smiles': 'O[C@@H]([C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)N[C@@H](CC(C)C)C(O)=O)C',
                                     'name': 'Glu-Thr-Leu',
                                     'reason': 'No benzene ring found'},
                                 {   'smiles': 'CCOc1ccc(NC(=O)C(C)O)cc1',
                                     'name': 'p-Lactophenetide',
                                     'reason': 'No benzene ring with exactly 4 '
                                               'chlorines found'},
                                 {   'smiles': 'O=C1NCC=2C1=C3C(N([C@@H]4O[C@H]([C@@H](O)[C@H]([C@H]4OC)O)C)C5=C3C=CC=C5)=C6NC7=C(C26)C=CC=C7',
                                     'name': "3'-epi-5'-methoxy-K252d",
                                     'reason': 'No benzene ring with exactly 4 '
                                               'chlorines found'},
                                 {   'smiles': 'O=C(OC1C(O)C(OC(C1O)CO)OCC2OC(OCC(OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)C(O)C(C2O)O)CCCCCCCCCCCCCCC',
                                     'name': '[(2S)-2-hexadecanoyloxy-3-[(2S,3R,4S,5S,6R)-6-[[(2S,3R,4S,5S,6R)-4-hexadecanoyloxy-3,5-dihydroxy-6-(hydroxymethyl)tetrahydropyran-2-yl]oxymethyl]-3,4,5-trihydroxy-tetrahydropyran-2-yl]oxy-propyl] '
                                             '(9E,12E,15E)-octadeca-9,12,15-trienoate',
                                     'reason': 'No benzene ring found'},
                                 {   'smiles': 'O=C1C2=C(O)C=C(OC)C=C2C(=O)C3=C1[C@@H]([C@@H](O)[C@]([C@@H]3O)(O)C)C',
                                     'name': 'Altersolanol G',
                                     'reason': 'No benzene ring with exactly 4 '
                                               'chlorines found'},
                                 {   'smiles': '[H][C@]1(O[C@](O)(C[C@H](O)[C@H]1NC(=O)CO)C([O-])=O)[C@H](O)[C@H](O)CO',
                                     'name': 'N-glycoloyl-alpha-neuraminate',
                                     'reason': 'No benzene ring found'},
                                 {   'smiles': 'OC(C(O)C/C=C\\C/C=C\\C/C=C\\CC)C/C=C\\C/C=C\\CCC(O)=O',
                                     'name': '10,11-DiHDPE',
                                     'reason': 'No benzene ring found'},
                                 {   'smiles': '[Na+].[H][C@]12SCC(C)=C(N1C(=O)[C@H]2NC(=O)[C@H](N)c1ccccc1)C([O-])=O',
                                     'name': 'cephalexin sodium',
                                     'reason': 'No benzene ring with exactly 4 '
                                               'chlorines found'}],
    'sample_false_negatives': [],
    'attempt': 0,
    'reasoning': 'Let me help create a program to classify '
                 'tetrachlorobenzenes. Looking at the examples and definition, '
                 'we need to:\n'
                 '\n'
                 '1. Check for a benzene ring\n'
                 '2. Verify it has exactly 4 chlorine atoms attached to it\n'
                 '3. Handle cases where the tetrachlorobenzene is part of a '
                 'larger structure\n'
                 '\n'
                 "Here's the implementation:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 25,
    'num_false_positives': 25,
    'num_true_negatives': 142250,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.5,
    'recall': 1.0,
    'f1': 0.6666666666666666,
    'accuracy': 0.9998243148278285,
    'negative_predictive_value': 1.0}