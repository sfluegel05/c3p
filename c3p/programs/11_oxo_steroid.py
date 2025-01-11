"""
Classifies: CHEBI:47787 11-oxo steroid
"""
"""
Classifies: 11-oxo steroid
Definition: Any oxo steroid that has an oxo substituent at position 11
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_11_oxo_steroid(smiles: str):
    """
    Determines if a molecule is an 11-oxo steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (True if molecule is an 11-oxo steroid, reason for classification)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for basic steroid core (four fused rings)
    steroid_core = Chem.MolFromSmarts("[#6]1~[#6]~[#6]~[#6]2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~[#6]~4~[#6]~3~[#6]~2~[#6]~1")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # SMARTS pattern for 11-oxo group in steroid context
    # This pattern looks for the specific environment of the 11-position ketone
    # in the steroid ring system
    oxo_11_pattern = Chem.MolFromSmarts("[#6]1~[#6]~[#6]2~[#6]~[#6]~[#6]~3~[#6](=[O:1])~[#6]~[#6]~[#6]~[#6]~3~[#6]~2~[#6]~[#6]~1")
    
    # Find matches for the 11-oxo pattern
    matches = mol.GetSubstructMatches(oxo_11_pattern)
    if not matches:
        return False, "No ketone group at position 11"
    
    # Count carbons to verify it's in the typical steroid range
    carbon_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6])
    if carbon_count < 19 or carbon_count > 30:
        return False, f"Carbon count ({carbon_count}) outside typical steroid range (19-30)"
        
    # Additional check for reasonable molecular weight
    mol_wt = Chem.Descriptors.ExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 500:
        return False, f"Molecular weight ({mol_wt:.1f}) outside typical steroid range (250-500)"
        
    # Count rings to ensure we have the right ring system
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings for steroid structure"
        
    return True, "Molecule contains steroid core with ketone group at position 11"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:47787',
                          'name': '11-oxo steroid',
                          'definition': 'Any oxo steroid that has an oxo '
                                        'substituent at position 11.',
                          'parents': ['CHEBI:35789'],
                          'xrefs': ['KEGG:C01985'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': None,
    'sample_false_negatives': None,
    'attempt': 0,
    'reasoning': 'Let me help create a program to classify 11-oxo steroids. '
                 "Let's break this down:\n"
                 '\n'
                 '1. An 11-oxo steroid must have:\n'
                 '   - A steroid core structure (4 fused rings)\n'
                 '   - An oxo (=O) group at position 11\n'
                 '   - The basic steroid numbering system places position 11 '
                 'at a specific location in the ring system\n'
                 '\n'
                 '2. Key characteristics to check:\n'
                 '   - The presence of the steroid core (four fused rings)\n'
                 '   - The presence of a ketone group (C=O) at position 11\n'
                 '   - The overall carbon count should be consistent with '
                 'steroids (typically 19-27 carbons)\n'
                 '\n'
                 "Here's the program:",
    'success': False,
    'best': True,
    'error': "module 'rdkit.Chem' has no attribute 'Descriptors'",
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