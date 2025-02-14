"""
Classifies: CHEBI:33916 aldopentose
"""
"""
Classifies: CHEBI:28053 aldopentose
"""
from rdkit import Chem

def is_aldopentose(smiles: str):
    """
    Determines if a molecule is an aldopentose based on its SMILES string.
    An aldopentose is a monosaccharide with five carbon atoms and an aldehyde group at one end,
    which can exist in open-chain or cyclic forms (furanose or pyranose rings).
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is an aldopentose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count number of carbon atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons != 5:
        return False, f"Contains {num_carbons} carbon atoms, should be 5 for an aldopentose"
    
    # Check for open-chain form with aldehyde group at one end
    aldehyde_pattern = Chem.MolFromSmarts('[CX3H1](=O)[CH2][CH](O)[CH](O)CO')
    if mol.HasSubstructMatch(aldehyde_pattern):
        return True, "Contains open-chain aldopentose with terminal aldehyde group"
    
    # Check for furanose ring (5-membered ring with oxygen)
    furanose_pattern = Chem.MolFromSmarts('[C;!R]1OC([C;!R])([C;!R])C1[O]')
    if mol.HasSubstructMatch(furanose_pattern):
        # Ensure total carbons are 5 (excluding the ring oxygen)
        ring_atoms = mol.GetSubstructMatch(furanose_pattern)
        num_ring_carbons = sum(1 for idx in ring_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        if num_ring_carbons == 4:
            return True, "Contains furanose ring derived from aldopentose"
    
    # Check for pyranose ring (6-membered ring with oxygen)
    pyranose_pattern = Chem.MolFromSmarts('[C;!R]1OC([C;!R])([C;!R])[C;!R][C;!R]1[O]')
    if mol.HasSubstructMatch(pyranose_pattern):
        # Ensure total carbons are 5 (excluding the ring oxygen)
        ring_atoms = mol.GetSubstructMatch(pyranose_pattern)
        num_ring_carbons = sum(1 for idx in ring_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        if num_ring_carbons == 5:
            return True, "Contains pyranose ring derived from aldopentose"
    
    # Alternative open-chain pattern (less strict)
    general_aldopentose = Chem.MolFromSmarts('[O]=[C][C@H](O)[C@H](O)CO')
    if mol.HasSubstructMatch(general_aldopentose):
        return True, "Contains open-chain aldopentose structure"
    
    # If none of the patterns match, molecule is not an aldopentose
    return False, "Does not match aldopentose structures (must have five carbons and aldehyde group at end or form specific cyclic structures)"
    

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:28053',
                              'name': 'aldopentose',
                              'definition': 'A pentose with a (potential) aldehyde group at one end.',
                              'parents': ['CHEBI:16842', 'CHEBI:24869']},
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
        'attempt': 2,
        'success': True,
        'best': True,
        'error': '',
        'stdout': None}