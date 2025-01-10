"""
Classifies: CHEBI:2468 secondary alpha-hydroxy ketone
"""
"""
Classifies: secondary alpha-hydroxy ketone
"""

from rdkit import Chem

def is_secondary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is a secondary alpha-hydroxy ketone based on its SMILES string.
    A secondary alpha-hydroxy ketone contains a ketone group (C=O) and a hydroxyl group (OH)
    attached to the alpha carbon, which is secondary (attached to one hydrogen, one carbonyl carbon,
    one hydroxyl group, and one other organyl group).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES to molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate over all atoms in the molecule
    for atom in mol.GetAtoms():
        # Check for carbonyl carbon (C=O)
        if atom.GetAtomicNum() == 6:
            is_carbonyl = False
            carbonyl_oxygen = None
            for neighbor in atom.GetNeighbors():
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                if neighbor.GetAtomicNum() == 8 and bond.GetBondType() == Chem.BondType.DOUBLE:
                    is_carbonyl = True
                    carbonyl_oxygen = neighbor
                    break
            if is_carbonyl:
                carbonyl_c = atom
                # Examine neighbors of the carbonyl carbon
                for alpha_c in carbonyl_c.GetNeighbors():
                    if alpha_c.GetIdx() == carbonyl_oxygen.GetIdx():
                        continue  # Skip the carbonyl oxygen
                    if alpha_c.GetAtomicNum() == 6:
                        # Check if alpha carbon has exactly one hydrogen
                        if alpha_c.GetTotalNumHs() == 1:
                            # Check if alpha carbon has hydroxyl group
                            has_oh = False
                            hydroxyl_oxygen = None
                            for nbr in alpha_c.GetNeighbors():
                                bond = mol.GetBondBetweenAtoms(alpha_c.GetIdx(), nbr.GetIdx())
                                if nbr.GetAtomicNum() == 8 and bond.GetBondType() == Chem.BondType.SINGLE:
                                    # Check if oxygen is bonded to hydrogen (hydroxyl group)
                                    if nbr.GetTotalNumHs() == 1:
                                        has_oh = True
                                        hydroxyl_oxygen = nbr
                                        break
                            if has_oh:
                                # Check if alpha carbon is connected to another organyl group
                                organyl_groups = [nbr for nbr in alpha_c.GetNeighbors() if nbr.GetIdx() not in [carbonyl_c.GetIdx(), hydroxyl_oxygen.GetIdx()]]
                                if len(organyl_groups) == 1:
                                    return True, "Molecule matches secondary alpha-hydroxy ketone pattern"
    # If no matching functional group is found
    return False, "Molecule does not match secondary alpha-hydroxy ketone pattern"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'secondary alpha-hydroxy ketone',
        'definition': 'An alpha-hydroxy ketone in which the carbonyl group and the hydroxy group '
                      'are linked by a carbon bearing one hydrogen and one organyl group. Secondary '
                      'alpha-hydroxy ketones are also known as acyloins, and are formally derived '
                      'from reductive coupling of two carboxylic acid groups.',
        'parents': []
    },
    'config': {
        'llm_model_name': 'your_model_name_here',
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
    'success': None,
    'best': None,
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