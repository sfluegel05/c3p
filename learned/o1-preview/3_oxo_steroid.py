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
    
    # Define the steroid backbone pattern (cyclopenta[a]phenanthrene ring system)
    steroid_pattern = Chem.MolFromSmarts('C1CCC2C3C(C1)CCC3CCC2')
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone detected"
    
    # Check for oxo group (ketone) at position 3
    # Define pattern for ketone attached to carbon at position 3 of steroid backbone
    # Since matching exact positions is complex, we will look for C=O attached to ring system
    ketone_pattern = Chem.MolFromSmarts('C(=O)[C]')  # Ketone group
    matches = mol.GetSubstructMatches(ketone_pattern)
    if not matches:
        return False, "No ketone group found"

    # Check if ketone is at position 3
    # Map the atoms in the steroid backbone to identify position 3
    # This requires a more detailed pattern matching
    steroid_core = Chem.MolFromSmarts("""
        [#6]1
        [#6][#6][#6]2
        [#6]3
        [#6]([#6]1)
        [#6][#6][#6]3
        [#6][#6][#6]2
    """)
    matches = mol.GetSubstructMatch(steroid_core)
    if not matches:
        return False, "Could not map steroid core for position numbering"
    else:
        # Assuming atom indices correspond to positions, get atom at position 3
        # Note: RDKit atom indices start from 0
        try:
            pos3_atom_idx = matches[2]  # Index of atom at position 3 in the steroid core
            pos3_atom = mol.GetAtomWithIdx(pos3_atom_idx)
            # Check if the atom at position 3 is connected to a ketone group
            ketone_at_pos3 = False
            for neighbor in pos3_atom.GetNeighbors():
                bond = mol.GetBondBetweenAtoms(pos3_atom_idx, neighbor.GetIdx())
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and neighbor.GetAtomicNum() == 8:
                    ketone_at_pos3 = True
                    break
            if not ketone_at_pos3:
                return False, "No ketone group at position 3"
        except IndexError:
            return False, "Error accessing atom at position 3"

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