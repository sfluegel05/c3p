"""
Classifies: CHEBI:4194 D-hexose
"""
"""
Classifies: CHEBI:15903 D-hexose
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_D_hexose(smiles: str):
    """
    Determines if a molecule is a D-hexose based on its SMILES string.
    A D-hexose is a hexose that has D-configuration at position 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a D-hexose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule has a ring structure
    if not mol.GetRingInfo().NumRings():
        return False, "No ring structure found"

    # Check for the presence of multiple hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[C][O]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 4:
        return False, "Not enough hydroxyl groups for a hexose"

    # Find the carbon at position 5 (second-to-last carbon in the ring)
    # We assume that the ring is a 6-membered ring (typical for hexoses)
    ring_atoms = mol.GetRingInfo().AtomRings()[0]
    if len(ring_atoms) != 6:
        return False, "Not a 6-membered ring"

    # Position 5 is the second-to-last carbon in the ring
    carbon_at_position_5 = mol.GetAtomWithIdx(ring_atoms[-2])

    # Check if the carbon at position 5 has a hydroxyl group
    has_hydroxyl = False
    for neighbor in carbon_at_position_5.GetNeighbors():
        if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 1:
            has_hydroxyl = True
            break

    if not has_hydroxyl:
        return False, "No hydroxyl group at position 5"

    # Check the stereochemistry of the carbon at position 5
    # In RDKit, the chirality is represented by the ChiralTag
    if carbon_at_position_5.GetChiralTag() != Chem.ChiralType.CHI_TETRAHEDRAL_CW:
        return False, "Carbon at position 5 does not have D-configuration"

    return True, "Hexose with D-configuration at position 5"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:15903',
                          'name': 'D-hexose',
                          'definition': 'A hexose that has D-configuration at position 5.',
                          'parents': ['CHEBI:15903', 'CHEBI:15903']},
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
    'stdout': None,
    'num_true_positives': 150,
    'num_false_positives': 4,
    'num_true_negatives': 182407,
    'num_false_negatives': 23,
    'num_negatives': None,
    'precision': 0.974025974025974,
    'recall': 0.8670520231213873,
    'f1': 0.9174311926605504,
    'accuracy': 0.9998521228585199}