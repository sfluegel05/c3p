"""
Classifies: CHEBI:35910 2-oxo monocarboxylic acid
"""
"""
Classifies: 2-oxo monocarboxylic acid
"""

from rdkit import Chem

def is_2_oxo_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid based on its SMILES string.
    A 2-oxo monocarboxylic acid is a molecule that contains at least one monocarboxylic acid group
    with a ketone (oxo) group at the alpha carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-oxo monocarboxylic acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1]")  # Carboxylic acid group
    ketone_pattern = Chem.MolFromSmarts("C=O")  # Ketone group

    # Find all carboxylic acid groups
    carboxy_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxy_matches:
        return False, "No carboxylic acid groups found"

    # Iterate over all carboxylic acid groups
    for carboxy_match in carboxy_matches:
        carboxyl_c_idx = carboxy_match[0]  # Index of the carbon in the carboxylic acid
        carboxyl_c_atom = mol.GetAtomWithIdx(carboxyl_c_idx)

        # Find alpha carbons (neighboring carbons to carboxyl carbon that are carbons)
        alpha_c_atoms = [nbr for nbr in carboxyl_c_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]

        # Check each alpha carbon for a ketone group
        for alpha_c_atom in alpha_c_atoms:
            alpha_c_idx = alpha_c_atom.GetIdx()

            # Check for ketone group attached to alpha carbon
            for bond in alpha_c_atom.GetBonds():
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    other_atom = bond.GetOtherAtom(alpha_c_atom)
                    if other_atom.GetAtomicNum() == 8:  # Oxygen atom
                        # Verify that this oxygen is part of a ketone group
                        # Ensure it's not part of the carboxylic acid group itself
                        if other_atom.GetIdx() != carboxy_match[1] and other_atom.GetIdx() != carboxy_match[2]:
                            return True, "Contains monocarboxylic acid group with oxo group at alpha carbon (2-oxo)"

            # Alternatively, use SMARTS to find ketone attached to alpha carbon
            ketone_matches = alpha_c_atom.Match(ketone_pattern)
            if ketone_matches:
                return True, "Contains monocarboxylic acid group with oxo group at alpha carbon (2-oxo)"

    # If no alpha carbon with ketone group is found
    return False, "No 2-oxo monocarboxylic acid substructure found"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': '2-oxo monocarboxylic acid',
        'definition': 'Any monocarboxylic acid having a 2-oxo substituent.'
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
    'attempt': 3,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
}