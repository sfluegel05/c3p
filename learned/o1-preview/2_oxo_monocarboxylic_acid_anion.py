"""
Classifies: CHEBI:35179 2-oxo monocarboxylic acid anion
"""
"""
Classifies: CHEBI:25450 2-oxo monocarboxylic acid anion
"""
from rdkit import Chem

def is_2_oxo_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid anion based on its SMILES string.
    A 2-oxo monocarboxylic acid anion is an oxo monocarboxylic acid anion in which the oxo group is located at the 2-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-oxo monocarboxylic acid anion, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylate groups and count them
    # Carboxylate group: -C(=O)[O-]
    carboxylate_pattern = Chem.MolFromSmarts("[C;H0](=O)[O-]")
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    num_carboxylate_groups = len(carboxylate_matches)
    if num_carboxylate_groups != 1:
        return False, f"Molecule has {num_carboxylate_groups} carboxylate groups, expected exactly one"

    # Identify the carboxylate carbon atom index
    carboxylate_C_idx = carboxylate_matches[0][0]
    carboxylate_C = mol.GetAtomWithIdx(carboxylate_C_idx)

    # Find alpha carbons (adjacent carbons to carboxylate carbon)
    alpha_carbons = []
    for neighbor in carboxylate_C.GetNeighbors():
        if neighbor.GetAtomicNum() == 6:  # Carbon atom
            alpha_carbons.append(neighbor)

    if not alpha_carbons:
        return False, "No alpha carbon adjacent to carboxylate group"

    # Check each alpha carbon for an oxo group at the 2-position
    for alpha_C in alpha_carbons:
        alpha_C_idx = alpha_C.GetIdx()
        # Look for oxo group attached to alpha carbon
        oxo_pattern = Chem.MolFromSmarts("[C;H0;$([C](=O))]")
        oxo_match = mol.GetSubstructMatches(oxo_pattern)
        for match in oxo_match:
            if alpha_C_idx == match[0]:
                return True, "Molecule has a carboxylate group with an oxo group at the 2-position"
        # Also check for ketal structures or acetyl groups at alpha carbon
        # Look for alpha carbon connected to a carbonyl carbon (e.g., acetyl group)
        carbonyl_pattern = Chem.MolFromSmarts("[C;H0]-[C](=O)-[C]")
        carbonyl_match = mol.GetSubstructMatches(carbonyl_pattern)
        for match in carbonyl_match:
            if alpha_C_idx == match[0]:
                return True, "Molecule has a carbonyl substituent at the 2-position"

    return False, "No oxo group found at the 2-position relative to the carboxylate group"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:25450',
        'name': '2-oxo monocarboxylic acid anion',
        'definition': 'An oxo monocarboxylic acid anion in which the oxo group is located at the 2-position.',
        'parents': ['CHEBI:35756', 'CHEBI:37549']
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
    }
}