"""
Classifies: CHEBI:47788 3-oxo steroid
"""
"""
Classifies: CHEBI:13653 3-oxo steroid
"""
from rdkit import Chem

def is_3_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo steroid based on its SMILES string.
    A 3-oxo steroid is any oxo steroid where an oxo substituent (C=O) is located at position 3.
    The function detects the steroid backbone and checks for a ketone at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the steroid nucleus with atom mapping to positions
    steroid_core_smarts = """
    [#6]1([#6H])-[#6]2-[#6]3([#6H])-[#6H]-[#6]-[#6]-4-[#6H]-[#6]-[#6]-[#6]-[#6]-[#6]-4-[#6]-3-[#6H]-[#6]-2-[#6]-1
    """
    steroid_core = Chem.MolFromSmarts(steroid_core_smarts)
    if steroid_core is None:
        return False, "Invalid steroid core SMARTS pattern"
    
    # Find the steroid core in the molecule
    match = mol.GetSubstructMatch(steroid_core)
    if not match:
        return False, "No steroid backbone detected"
    
    # Map positions based on the steroid core
    # Standard steroid numbering has positions 1 to 17
    # For accurate mapping, we should use a SMARTS pattern with atom mapping numbers
    steroid_numbering_smarts = """
    [$([#6]1)](-[#6H])-[#6]2(-[#6H])-[#6](=[#6])-[$([#6H])]-[#6]3=[#6]-[#6]-[#6]-[#6]-[#6]-3-[#6]-2-[#6]-1
    """
    steroid_numbering = Chem.MolFromSmarts(steroid_numbering_smarts)
    steroid_matches = mol.GetSubstructMatches(steroid_numbering)
    
    if not steroid_matches:
        return False, "Steroid core not matched with numbering"

    # Assuming the first match is the steroid nucleus
    atom_indices = steroid_matches[0]
    
    # Position 3 is the third atom in the SMARTS pattern (index 2 in atom_indices)
    pos3_atom_idx = atom_indices[2]
    pos3_atom = mol.GetAtomWithIdx(pos3_atom_idx)
    if pos3_atom.GetAtomicNum() != 6:
        return False, "Atom at position 3 is not carbon"
    
    # Check for a ketone (C=O) at position 3
    ketone_pattern = Chem.MolFromSmarts("[#6]=O")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    ketone_atom_indices = [match[0] for match in ketone_matches]
    if pos3_atom_idx not in ketone_atom_indices:
        return False, "No ketone group at position 3"
    
    return True, "Molecule is a 3-oxo steroid with ketone at position 3"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:13653',
        'name': '3-oxo steroid',
        'definition': 'Any oxo steroid where an oxo substituent is located at position 3.',
        'parents': ['CHEBI:36858', 'CHEBI:36686']
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': 50,
        'max_negative_in_prompt': 20,
        'max_instances_in_prompt': 100,
        'test_proportion': 0.1
    },
    'message': None,
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None
}