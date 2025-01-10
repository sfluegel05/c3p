"""
Classifies: CHEBI:33913 corrinoid
"""
"""
Classifies: CHEBI:33913 corrinoid
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_corrinoid(smiles: str):
    """
    Determines if a molecule is a corrinoid based on its SMILES string.
    A corrinoid is a derivative of the corrin nucleus, which contains four reduced or partly reduced
    pyrrole rings joined in a macrocycle by three =C- groups and one direct carbon-carbon bond linking
    alpha positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a corrinoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Simplified corrin nucleus SMARTS pattern
    # The corrin nucleus is a macrocycle with four pyrrole rings connected via three methine bridges (=C-)
    # and one direct carbon-carbon bond between two pyrrole alpha positions.

    # Due to the complexity, we create a simplified representation of the corrin ring system.
    corrin_smarts = """
    C1=C[C@H]2C=C[C@@H]3C=C[C@H]4C=CC[C@@H]1N2[Co]N3C4
    """

    # Convert SMARTS to mol
    corrin_pattern = Chem.MolFromSmarts(corrin_smarts)
    if corrin_pattern is None:
        return False, "Error in SMARTS pattern for corrin nucleus"

    # Check for substructure match
    if mol.HasSubstructMatch(corrin_pattern):
        return True, "Molecule contains the corrin nucleus"
    else:
        return False, "Molecule does not contain the corrin nucleus"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:33913',
        'name': 'corrinoid',
        'definition': 'A derivative of the corrin nucleus, which contains four reduced or partly reduced pyrrole rings joined in a macrocycle by three =C- groups and one direct carbon-carbon bond linking alpha positions.',
        'parents': []
    },
    'config': {
        # Configuration parameters can be included here if necessary
    },
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    # Metrics placeholders
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