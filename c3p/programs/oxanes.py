"""
Classifies: CHEBI:46942 oxanes
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_oxanes(smiles: str):
    """
    Determines if a molecule is an oxane or its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oxane or its substituted derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()

    # Find all 6-membered rings
    six_membered_rings = [ring for ring in rings.AtomRings() if len(ring) == 6]

    if not six_membered_rings:
        return False, "No 6-membered rings found"

    # Check for presence of oxygen in the ring
    for ring in six_membered_rings:
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if any(atom.GetSymbol() == 'O' for atom in atoms):
            return True, "Contains an oxane ring"

    return False, "No oxane rings found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:46942',
                          'name': 'oxanes',
                          'definition': 'Any organic heteromonocyclic '
                                        'compoundthat is oxane or its '
                                        'substituted derivatives.',
                          'parents': ['CHEBI:25693', 'CHEBI:38104']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 25,
    'num_false_positives': 2,
    'num_true_negatives': 18,
    'num_false_negatives': 2,
    'precision': 0.9259259259259259,
    'recall': 0.9259259259259259,
    'f1': 0.9259259259259259,
    'accuracy': None}