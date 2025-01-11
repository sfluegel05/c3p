"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
"""
Classifies: CHEBI:26873 sesquiterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.
    Sesquiterpenoids are terpenoids derived from sesquiterpenes (C15 skeleton),
    which may be rearranged or modified by the removal of methyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesquiterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 13 or c_count > 15:
        return False, f"Number of carbons ({c_count}) not in range 13-15 for sesquiterpenoids"

    # Count isoprene units (C5 units)
    # Define SMARTS pattern for isoprene unit
    isoprene_pattern = Chem.MolFromSmarts('C(=C)C=C')  # Simplified isoprene pattern
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 1:
        return False, "No isoprene units found"

    # Check for terpenoid functional groups
    # Common functional groups: hydroxyl, ketone, aldehyde, ester, ether
    functional_groups = ['[OX2H]', '[CX3]=[OX1]', '[CX3H1](=O)', '[CX3](=O)[OX2H1]', '[OX2][CX4]']
    has_functional_group = False
    for fg in functional_groups:
        fg_pattern = Chem.MolFromSmarts(fg)
        if mol.HasSubstructMatch(fg_pattern):
            has_functional_group = True
            break
    if not has_functional_group:
        return False, "No common terpenoid functional groups found"

    # Optionally, check for ring structures
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 0:
        # Sesquiterpenoids can be acyclic, so we don't fail here
        ring_message = "Molecule is acyclic"
    else:
        ring_message = f"Molecule has {ring_info.NumRings()} ring(s)"

    return True, f"Classified as sesquiterpenoid: {c_count} carbons, contains isoprene units, {ring_message}"

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26873',
                              'name': 'sesquiterpenoid',
                              'definition': 'Any terpenoid derived from a sesquiterpene. The term includes compounds in which the C15 skeleton of the parent sesquiterpene has been rearranged or modified by the removal of one or more skeletal atoms (generally methyl groups).',
                              'parents': ['CHEBI:26870', 'CHEBI:48190']},
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