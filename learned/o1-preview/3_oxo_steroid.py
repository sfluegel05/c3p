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
    steroid_core = Chem.MolFromSmarts("""
    [#6]12[#6]3[#6](~[#6]~[#6]~[#6]~1)~[#6]~[#6]~[#6]~[#6]~2~[#6]~[#6]~3
    """)
    if steroid_core is None:
        return False, "Invalid steroid core SMARTS pattern"

    # Find the steroid core in the molecule
    match = mol.GetSubstructMatch(steroid_core)
    if not match:
        return False, "No steroid backbone detected"
    
    # Map positions based on the steroid core
    # Assuming atom indices in match correspond to positions 1 to 17
    # This is a simplification; actual mapping may require more precise definition
    pos3_atom_idx = match[2]  # Index of the atom at position 3
    
    # Check if the atom at position 3 is a carbon
    pos3_atom = mol.GetAtomWithIdx(pos3_atom_idx)
    if pos3_atom.GetAtomicNum() != 6:
        return False, "Atom at position 3 is not carbon"
    
    # Check for a ketone (C=O) at position 3
    ketone_pattern = Chem.MolFromSmarts("[#6](=O)")

    # Check if the carbon at position 3 has a double bond to oxygen
    has_ketone = False
    for bond in pos3_atom.GetBonds():
        neighbor = bond.GetOtherAtom(pos3_atom)
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and neighbor.GetAtomicNum() == 8:
            has_ketone = True
            break
    if not has_ketone:
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
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None
}