"""
Classifies: CHEBI:47788 3-oxo steroid
"""
"""
Classifies: CHEBI:13653 3-oxo steroid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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
    
    # Define a general steroid backbone pattern (cyclopentanoperhydrophenanthrene)
    steroid_pattern = Chem.MolFromSmarts("""
        [#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1                          # Aromatic A ring
        -[#6]2-[#6]-[#6]-[#6]-3-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-3   # Rings B and C
        -[#6]-2                                                  # Ring D
    """)
    # Alternatively, use a more general pattern for the steroid nucleus
    steroid_pattern = Chem.MolFromSmarts("""
        [#6]1[#6][#6][#6]2[#6]1[#6][#6][#6]3[#6][#6][#6][#6][#6]23  # Steroid core
    """)
    matches = mol.GetSubstructMatches(steroid_pattern)
    if not matches:
        return False, "No steroid backbone detected"

    # Map the atoms in the steroid core to standard numbering
    # Define the steroid core with atom mapping numbers
    steroid_core = Chem.MolFromSmarts("""
        [#6]1([#6])[#6][#6]2[#6](1)[#6][#6]3[#6](2)[#6][#6][#6][#6]3
    """)
    match = mol.GetSubstructMatch(steroid_core)
    if not match:
        return False, "Could not map steroid core for position numbering"
    
    # Create a mapping from atom map numbers to atom indices in the molecule
    atom_map = {}
    for pattern_atom_idx, mol_atom_idx in enumerate(match):
        atom_map[pattern_atom_idx + 1] = mol_atom_idx  # Atom map numbers start from 1
    
    # Identify atom at position 3 in the steroid core
    pos3_atom_idx = atom_map.get(3)
    if pos3_atom_idx is None:
        return False, "Could not identify atom at position 3"

    # Check if the atom at position 3 is connected to a ketone group (C=O)
    pos3_atom = mol.GetAtomWithIdx(pos3_atom_idx)
    ketone_at_pos3 = False
    for neighbor in pos3_atom.GetNeighbors():
        bond = mol.GetBondBetweenAtoms(pos3_atom_idx, neighbor.GetIdx())
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and neighbor.GetAtomicNum() == 8:
            ketone_at_pos3 = True
            break
    if not ketone_at_pos3:
        return False, "No ketone group at position 3"

    return True, "Molecule is a 3-oxo steroid with ketone at position 3"

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:13653',
                              'name': '3-oxo steroid',
                              'definition': 'Any oxo steroid where an oxo '
                                            'substituent is located at position '
                                            '3.',
                              'parents': ['CHEBI:36858', 'CHEBI:36686']},
        'config': {   'llm_model_name': 'lbl/claude-sonnet',
                      'f1_threshold': 0.8,
                      'max_attempts': 5,
                      'max_positive_instances': None,
                      'max_positive_to_test': None,
                      'max_negative_to_test': None,
                      'max_positive_in_prompt': 50,
                      'max_negative_in_prompt': 20,
                      'max_instances_in_prompt': 100,
                      'test_proportion': 0.1},
        'message': None,
        'attempt': 0,
        'success': True,
        'best': True,
        'error': '',
        'stdout': None}