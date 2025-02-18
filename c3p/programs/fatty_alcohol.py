"""
Classifies: CHEBI:24026 fatty alcohol
"""
"""
Classifies: CHEBI:30843 fatty alcohol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a fatty alcohol based on its SMILES string.
    A fatty alcohol is an aliphatic alcohol with a carbon chain length of 3 to 27 atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for at least one alcohol group (-OH)
    alcohol_pattern = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(alcohol_pattern):
        return False, "No alcohol group found"

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 3 or c_count > 27:
        return False, f"Carbon chain length {c_count} is outside the range of 3 to 27"

    # Ensure the molecule is aliphatic (no aromatic rings)
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    if aromatic_atoms > 0:
        return False, "Molecule contains aromatic rings, not aliphatic"

    # Check if the molecule is primarily a hydrocarbon chain with an alcohol group
    # This is a heuristic and may need refinement
    h_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count != 1:
        return False, f"Expected exactly 1 oxygen (alcohol group), found {o_count}"

    # Check for unsaturation (double or triple bonds)
    unsaturation = rdMolDescriptors.CalcNumUnsaturatedBonds(mol)
    if unsaturation > 0:
        # Allow unsaturated fatty alcohols
        pass

    # Check for branching
    # This is optional since fatty alcohols can be branched or unbranched

    return True, f"Aliphatic alcohol with {c_count} carbon atoms and 1 alcohol group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:30843',
                          'name': 'fatty alcohol',
                          'definition': 'An aliphatic alcohol consisting of a '
                                        'chain of 3 to greater than 27 carbon '
                                        'atoms. Fatty alcohols may be saturated '
                                        'or unsaturated and may be branched or '
                                        'unbranched.',
                          'parents': ['CHEBI:30843', 'CHEBI:25805']},
    'config': {   'llm_model_name': 'lbl/claude-sonnet',
                  'f1_threshold': 0.8,
                  'max_attempts': 5,
                  'max_positive_instances': None,
                  'max_positive_to_test': None,
                  'max_negative_to_test': None,
                  'max_positive_in_prompt': 50,
                  'max_negative_in_prompt': 20,
                  'max_instances_in_prompt': 100,
                  'test_proportion': 0.1},
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 150,
    'num_false_positives': 4,
    'num_true_negatives': 182407,
    'num_false_negatives': 23,
    'num_negatives': None,
    'precision': 0.974025974025974,
    'recall': 0.8670520231213873,
    'f1': 0.9174311926605504,
    'accuracy': 0.9998521228585199}