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
    A fatty alcohol is an aliphatic alcohol with a chain of 3 to greater than 27 carbon atoms.
    The chain can be saturated or unsaturated and may be branched or unbranched.

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

    # Check for the presence of at least one hydroxyl group (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group found"

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 3:
        return False, f"Too few carbon atoms ({c_count}), need at least 3"
    if c_count > 27:
        return True, f"Chain length ({c_count}) is greater than 27 carbon atoms"

    # Check if the molecule is aliphatic (no aromatic rings)
    aromatic_atoms = [atom for atom in mol.GetAtoms() if atom.GetIsAromatic()]
    if aromatic_atoms:
        return False, "Molecule contains aromatic rings, not aliphatic"

    # Check if the molecule is a chain (no cycles)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings, not a chain"

    # Check if the molecule is branched or unbranched
    # This is more complex, but we can check the number of carbon atoms with more than 2 connections
    branched_carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetDegree() > 2]
    if branched_carbons:
        return True, f"Branched chain with {c_count} carbon atoms"
    else:
        return True, f"Unbranched chain with {c_count} carbon atoms"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:30843',
                          'name': 'fatty alcohol',
                          'definition': 'An aliphatic alcohol consisting of a '
                                        'chain of 3 to greater than 27 carbon '
                                        'atoms. Fatty alcohols may be saturated '
                                        'or unsaturated and may be branched or '
                                        'unbranched.',
                          'parents': ['CHEBI:30843', 'CHEBI:30843']},
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