"""
Classifies: CHEBI:27933 beta-lactam antibiotic
"""
"""
Classifies: beta-lactam antibiotic
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_beta_lactam_antibiotic(smiles: str):
    """
    Determines if a molecule is a beta-lactam antibiotic based on its SMILES string.
    A beta-lactam antibiotic is an organonitrogen heterocyclic antibiotic that contains a beta-lactam ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-lactam antibiotic, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for beta-lactam ring
    # Beta-lactam ring: 4-membered cyclic amide (lactam) containing a nitrogen atom and a carbonyl group
    beta_lactam_smarts = '[R4;$([N&R][C&R][C&R][C&R](=O))]'
    beta_lactam_pattern = Chem.MolFromSmarts(beta_lactam_smarts)
    if beta_lactam_pattern is None:
        return False, "Error in beta-lactam SMARTS pattern"

    # Search for beta-lactam ring
    if not mol.HasSubstructMatch(beta_lactam_pattern):
        return False, "No beta-lactam ring found"

    return True, "Contains a beta-lactam ring"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'beta-lactam antibiotic',
        'definition': 'An organonitrogen heterocyclic antibiotic that contains a beta-lactam ring.',
        'parents': []
    },
    'config': {},
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
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