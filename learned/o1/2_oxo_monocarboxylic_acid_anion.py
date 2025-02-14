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

    # Check for carboxylate groups (-C(=O)[O-])
    carboxylate_pattern = Chem.MolFromSmarts("[C](=O)[O-]")
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    if not carboxylate_matches:
        return False, "No carboxylate groups found"

    # For each carboxylate group, check for an oxo group at the alpha carbon
    for match in carboxylate_matches:
        carboxylate_C_idx = match[0]
        carboxylate_C = mol.GetAtomWithIdx(carboxylate_C_idx)

        # Find alpha carbons (adjacent carbons)
        alpha_carbons = [nbr for nbr in carboxylate_C.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if not alpha_carbons:
            continue  # No alpha carbon, move to next carboxylate group

        for alpha_C in alpha_carbons:
            # Check if alpha_C has a double bond to oxygen (oxo group)
            oxo_found = False
            for bond in alpha_C.GetBonds():
                neighbor = bond.GetOtherAtom(alpha_C)
                if neighbor.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    oxo_found = True
                    break
            if oxo_found:
                return True, "Molecule is a 2-oxo monocarboxylic acid anion"

    return False, "No alpha carbon with oxo group found adjacent to any carboxylate group"

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
    },
    'message': None,
    'attempt': 1,
    'success': False,
    'best': False,
    'error': '',
    'stdout': None,
    'num_true_positives': None,
    'num_false_positives': None,
    'num_true_negatives': None,
    'num_false_negatives': None,
    'num_negatives': None,
    'precision': None,
    'recall': None,
    'f1': None,
    'accuracy': None
}