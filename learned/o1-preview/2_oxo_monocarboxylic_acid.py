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
    A 2-oxo monocarboxylic acid is a monocarboxylic acid with a ketone group at the 2-position (alpha carbon).

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

    # Check for exactly one carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1]")
    carboxy_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxy_matches) != 1:
        return False, f"Found {len(carboxy_matches)} carboxylic acid groups, need exactly 1"

    # Get the carboxyl carbon atom index
    carboxyl_c_idx = carboxy_matches[0][0]
    carboxyl_c_atom = mol.GetAtomWithIdx(carboxyl_c_idx)

    # Find alpha carbon (carbon neighbor of carboxyl carbon)
    alpha_c_atoms = [atom for atom in carboxyl_c_atom.GetNeighbors() if atom.GetAtomicNum() == 6]
    if len(alpha_c_atoms) == 0:
        return False, "No alpha carbon connected to carboxyl carbon"

    # Check if alpha carbon has a ketone group (=O)
    ketone_found = False
    for alpha_c_atom in alpha_c_atoms:
        # Check for double bond to oxygen on alpha carbon
        for bond in alpha_c_atom.GetBonds():
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                other_atom = bond.GetOtherAtom(alpha_c_atom)
                if other_atom.GetAtomicNum() == 8:  # Oxygen atom
                    ketone_found = True
                    break
        if ketone_found:
            break

    if not ketone_found:
        return False, "No ketone group at alpha carbon"

    return True, "Contains monocarboxylic acid group with ketone at alpha carbon (2-oxo)"

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
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
}