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

    # Identify carboxylate groups (-C(=O)[O-]) or carboxylic acids (-C(=O)O)
    carboxylate_pattern = Chem.MolFromSmarts("[CX3](=O)[OX1-]")
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX1H]")
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    carboxy_matches = carboxylate_matches + carboxylic_acid_matches

    num_carboxy_groups = len(carboxy_matches)
    if num_carboxy_groups == 0:
        return False, "No carboxylate or carboxylic acid groups found"

    # For each carboxylate/carboxylic acid group, check for oxo group at 2-position
    for match in carboxy_matches:
        carboxy_C_idx = match[0]  # Index of the carbonyl carbon
        carboxy_C = mol.GetAtomWithIdx(carboxy_C_idx)

        # Find alpha carbons (adjacent carbons to carboxylate carbon)
        alpha_carbons = []
        for neighbor in carboxy_C.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon atom
                alpha_carbons.append(neighbor)

        if not alpha_carbons:
            continue  # No alpha carbon adjacent to this carboxylate group

        # Check each alpha carbon for an oxo group at the 2-position
        for alpha_C in alpha_carbons:
            alpha_C_idx = alpha_C.GetIdx()

            # Check if alpha carbon is double-bonded to an oxygen (oxo group)
            has_oxo = False
            for neighbor in alpha_C.GetNeighbors():
                if neighbor.GetAtomicNum() == 8:  # Oxygen atom
                    bond = mol.GetBondBetweenAtoms(alpha_C_idx, neighbor.GetIdx())
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        has_oxo = True
                        break

            if has_oxo:
                return True, "Molecule has a carboxylate group with an oxo group at the 2-position"

    return False, "No oxo group found at the 2-position relative to any carboxylate group"

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