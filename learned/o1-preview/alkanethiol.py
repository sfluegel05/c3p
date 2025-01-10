"""
Classifies: CHEBI:47908 alkanethiol
"""
"""
Classifies: CHEBI:29255 alkanethiol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alkanethiol(smiles: str):
    """
    Determines if a molecule is an alkanethiol based on its SMILES string.
    An alkanethiol is a compound in which a sulfanyl group (-SH) is attached to an alkyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanethiol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens to the molecule
    mol_with_H = Chem.AddHs(mol)

    # Flag to indicate whether an alkanethiol group is found
    has_alkanethiol = False

    # Iterate over atoms to find sulfur atoms
    for atom in mol_with_H.GetAtoms():
        if atom.GetAtomicNum() == 16:  # Sulfur atom
            neighbors = atom.GetNeighbors()
            if len(neighbors) != 2:
                continue  # Sulfur must have exactly 2 neighbors (H and C)
            num_hydrogens = 0
            num_carbons = 0
            for neighbor in neighbors:
                if neighbor.GetAtomicNum() == 1:
                    num_hydrogens += 1
                elif neighbor.GetAtomicNum() == 6:
                    num_carbons += 1
            if num_hydrogens == 1 and num_carbons == 1:
                has_alkanethiol = True
                break

    if has_alkanethiol:
        return True, "Contains an -SH group attached to an alkyl group"
    else:
        return False, "No -SH group attached to an alkyl group found"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:29255',
        'name': 'alkanethiol',
        'definition': 'An alkanethiol is a compound in which a sulfanyl group (-SH) is attached to an alkyl group.',
        'parents': ['CHEBI:24578']
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
    # Metrics would be filled in after evaluation
}