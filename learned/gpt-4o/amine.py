"""
Classifies: CHEBI:32952 amine
"""
from rdkit import Chem

def is_amine(smiles: str):
    """
    Determines if a molecule is an amine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate over atoms in the molecule
    for atom in mol.GetAtoms():
        # Check if the atom is Nitrogen
        if atom.GetAtomicNum() == 7:
            # Count the number of directly attached carbon atoms (hydrocarbyl groups)
            num_carbon = sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 6)
            # Count the number of directly attached hydrogen atoms
            num_hydrogen = sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 1)
            
            # Check if it matches primary, secondary, or tertiary amine
            if (num_carbon + num_hydrogen) == 3:
                if num_carbon > 0 and (3 - num_carbon) == num_hydrogen:
                    return True, f"Identified as an amine with {num_carbon} hydrocarbyl groups and {num_hydrogen} hydrogens attached to nitrogen"

    return False, "No amine group found"

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:32988',
                          'name': 'amine',
                          'definition': 'A compound formally derived from ammonia by replacing one, two or three hydrogen atoms by hydrocarbyl groups.'},
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
    'success': None,
    'best': None,
    'error': '',
    'stdout': None,
    'precision': None,
    'recall': None,
    'f1': None,
    'accuracy': None
}