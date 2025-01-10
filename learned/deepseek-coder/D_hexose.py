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

    # Check if the molecule has exactly 6 carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 6:
        return False, f"Not a hexose: found {c_count} carbon atoms, need exactly 6"

    # Check for the presence of a hydroxyl group at position 5 with D-configuration
    # We need to find a carbon with a hydroxyl group and check its chirality
    # The fifth carbon is typically the second-to-last carbon in the chain
    # We will look for a carbon with a hydroxyl group and check its stereochemistry

    # Define a pattern for a carbon with a hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts("[C][O]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    
    if not hydroxyl_matches:
        return False, "No hydroxyl groups found"

    # Check the stereochemistry of the hydroxyl group at position 5
    # We need to find the carbon at position 5 and check its chirality
    # The fifth carbon is typically the second-to-last carbon in the chain
    # We will look for a carbon with a hydroxyl group and check its stereochemistry

    # Get the atoms in the molecule
    atoms = mol.GetAtoms()
    
    # Find the carbon at position 5 (second-to-last carbon)
    # Assuming the molecule is linear, we can find the second-to-last carbon
    # with a hydroxyl group
    carbon_at_position_5 = None
    for atom in atoms:
        if atom.GetAtomicNum() == 6 and atom.GetDegree() == 4:
            # Check if this carbon has a hydroxyl group
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 1:
                    carbon_at_position_5 = atom
                    break
            if carbon_at_position_5:
                break

    if not carbon_at_position_5:
        return False, "No carbon at position 5 with a hydroxyl group found"

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