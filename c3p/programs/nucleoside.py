"""
Classifies: CHEBI:33838 nucleoside
"""
"""
Classifies: CHEBI:33838 nucleoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside(smiles: str):
    """
    Determines if a molecule is a nucleoside based on its SMILES string.
    A nucleoside consists of a nucleobase attached to ribose/deoxyribose via N-glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for furanose (5-membered sugar ring)
    furanose_pattern = Chem.MolFromSmarts("[CH2]-[CH1]-[CH1]-[CH1]-O1")
    if not mol.HasSubstructMatch(furanose_pattern):
        return False, "No ribose/deoxyribose sugar found"

    # Look for common nucleobase patterns
    base_patterns = [
        # Purine core
        "[nR1r6]1c[nR1r6]c2c1[nR1r6]c[nR1r6]c2",
        # Pyrimidine core 
        "[nR1r6]1c[nR1r6]cc(=O)[nR1r6]c1",
    ]
    
    found_base = False
    for pattern in base_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_base = True
            break
            
    if not found_base:
        return False, "No nucleobase pattern found"

    # Check for N-glycosidic bond (N-C1' bond between base and sugar)
    n_glycosidic = Chem.MolFromSmarts("[nR1r6]C1OC(CO)CC1")
    if not mol.HasSubstructMatch(n_glycosidic):
        return False, "No N-glycosidic bond found"

    # Basic size checks
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 800:
        return False, f"Molecular weight {mol_wt:.1f} outside typical range for nucleosides"

    # Count key atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 9:
        return False, "Too few carbons for nucleoside"
    if n_count < 2:
        return False, "Too few nitrogens for nucleoside"
    if o_count < 4:
        return False, "Too few oxygens for nucleoside"

    # Check for reasonable ring count (typically 2-3 rings)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 2:
        return False, "Too few rings for nucleoside structure"
    if ring_info.NumRings() > 4:
        return False, "Too many rings for typical nucleoside"

    return True, "Contains nucleobase attached to ribose/deoxyribose via N-glycosidic bond"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33838',
                          'name': 'nucleoside',
                          'definition': 'An N-glycosyl compound that has both '
                                        'a nucleobase, normally adenine, '
                                        'guanine, xanthine, thymine, cytosine '
                                        'or uracil, and either a ribose or '
                                        'deoxyribose as functional parents.',
                          'parents': [   'CHEBI:21731',
                                         'CHEBI:26912',
                                         'CHEBI:61120'],
                          'xrefs': ['KEGG:C00801', 'Wikipedia:Nucleoside'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': None,
    'sample_false_negatives': None,
    'attempt': 0,
    'reasoning': 'Let me help create a program to classify nucleosides. A '
                 'nucleoside consists of a nucleobase (like adenine, guanine, '
                 'cytosine, thymine, uracil) attached to a ribose or '
                 'deoxyribose sugar via an N-glycosidic bond.\n'
                 '\n'
                 "Here's my approach:\n"
                 '1. Check for valid SMILES and basic structure\n'
                 '2. Look for sugar moiety (ribose or deoxyribose)\n'
                 '3. Look for nucleobase patterns\n'
                 '4. Verify N-glycosidic bond connecting them\n'
                 '5. Additional checks for reasonable size and composition',
    'success': False,
    'best': True,
    'error': 'Python argument types in\n'
             '    Mol.HasSubstructMatch(Mol, NoneType)\n'
             'did not match C++ signature:\n'
             '    HasSubstructMatch(RDKit::ROMol self, RDKit::MolBundle query, '
             'RDKit::SubstructMatchParameters params=True)\n'
             '    HasSubstructMatch(RDKit::ROMol self, RDKit::ROMol query, '
             'RDKit::SubstructMatchParameters params)\n'
             '    HasSubstructMatch(RDKit::ROMol self, RDKit::MolBundle query, '
             'bool recursionPossible=True, bool useChirality=False, bool '
             'useQueryQueryMatches=False)\n'
             '    HasSubstructMatch(RDKit::ROMol self, RDKit::ROMol query, '
             'bool recursionPossible=True, bool useChirality=False, bool '
             'useQueryQueryMatches=False)',
    'stdout': '',
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None,
    'negative_predictive_value': 0.0}