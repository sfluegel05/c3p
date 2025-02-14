"""
Classifies: CHEBI:20156 3-oxo-Delta(1) steroid
"""
"""
Classifies: CHEBI:35164 3-oxo-Delta(1) steroid
"""
from rdkit import Chem

def is_3_oxo_Delta_1__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(1) steroid based on its SMILES string.
    A 3-oxo-Delta(1) steroid is any 3-oxo steroid that contains a double bond between positions 1 and 2.
    This function checks for:
      - Steroid backbone (three six-membered rings and one five-membered ring fused together)
      - Ketone group at position 3
      - Double bond between positions 1 and 2
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a 3-oxo-Delta(1) steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define steroid backbone pattern with atom mapping
    steroid_pattern = Chem.MolFromSmarts("""
    [#6]1([#6])[#6][#6]2[#6]([#6][#6]3[#6](=[#6][#6][#6]4[#6][#6][#6][#6][#6]4)[#6][#6]3)[#6][#6]12
    """)
    if steroid_pattern is None:
        return False, "Invalid steroid SMARTS pattern"
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Identify position 3 carbon and check for ketone group (C=O)
    ketone_at_3_pattern = Chem.MolFromSmarts("""
    [#6]1([#6])[#6][#6]2[#6]([#6][#6]3[#6](=[#6][#6][#6]4[#6][#6][#6][#6][#6]4)[#6][#6]3)[#6][#6]12
    :[#6]:[#6]:[#6]:[#6]:[#6]:[#6]-[#6](=O)-[#6]
    """)
    ketone_pattern = Chem.MolFromSmarts("[#6]-C(=O)")
    # Map the steroid backbone to get the atoms at positions
    match = mol.GetSubstructMatch(steroid_pattern)
    if not match:
        return False, "No steroid backbone found"
    # Positions 1 and 2 in the match
    atom1_idx = match[0]
    atom2_idx = match[1]
    atom3_idx = match[2]
    atom3 = mol.GetAtomWithIdx(atom3_idx)
    # Check if atom at position 3 has a ketone group
    ketone_found = False
    for bond in atom3.GetBonds():
        nbr = bond.GetOtherAtom(atom3)
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and nbr.GetAtomicNum() == 8:
            ketone_found = True
            break
    if not ketone_found:
        return False, "No ketone group at position 3"
    
    # Check for double bond between positions 1 and 2
    atom1 = mol.GetAtomWithIdx(atom1_idx)
    atom2 = mol.GetAtomWithIdx(atom2_idx)
    bond = mol.GetBondBetweenAtoms(atom1_idx, atom2_idx)
    if bond is None or bond.GetBondType() != Chem.rdchem.BondType.DOUBLE:
        return False, "No double bond between positions 1 and 2"
    
    return True, "Matches 3-oxo-Delta(1) steroid with correct backbone, ketone at position 3, and double bond between positions 1 and 2"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35164',
                              'name': '3-oxo-Delta(1) steroid',
                              'definition': 'Any 3-oxo steroid that contains a '
                                            'double bond between positions 1 and 2.'},
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
        'attempt': 1,
        'success': True,
        'best': True,
        'error': '',
        'stdout': None}