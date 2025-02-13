"""
Classifies: CHEBI:82744 aliphatic aldoxime
"""
"""
Classifies: CHEBI:36586 aliphatic aldoxime
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_aliphatic_aldoxime(smiles: str):
    """
    Determines if a molecule is an aliphatic aldoxime based on its SMILES string.
    An aliphatic aldoxime is derived from an aliphatic aldehyde and contains the aldoxime group (-CH=N-OH).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic aldoxime, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the aldoxime group pattern (-CH=N-OH)
    aldoxime_pattern = Chem.MolFromSmarts("[CX3H1]=[NX2][OX2H1]")
    if not mol.HasSubstructMatch(aldoxime_pattern):
        return False, "No aldoxime group (-CH=N-OH) found"

    # Check if the molecule is aliphatic (no aromatic rings or conjugated systems)
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Molecule contains aromatic rings, not aliphatic"

    # Check for conjugated systems (double bonds not part of the aldoxime group)
    conjugated_pattern = Chem.MolFromSmarts("[CX3]=[CX3]")
    conjugated_matches = mol.GetSubstructMatches(conjugated_pattern)
    if len(conjugated_matches) > 0:
        return False, "Molecule contains conjugated systems, not aliphatic"

    # Check for aliphatic chain length (at least one carbon in addition to the aldoxime group)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 2:
        return False, "Molecule does not have a sufficient aliphatic chain"

    return True, "Contains an aldoxime group (-CH=N-OH) attached to an aliphatic chain"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:36586',
        'name': 'aliphatic aldoxime',
        'definition': 'Any aldoxime derived from an aliphatic aldehyde.',
        'parents': ['CHEBI:36585', 'CHEBI:36584']
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
    'num_true_positives': 20,
    'num_false_positives': 0,
    'num_true_negatives': 100,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': 1.0
}