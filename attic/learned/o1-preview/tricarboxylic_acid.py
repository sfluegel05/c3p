"""
Classifies: CHEBI:27093 tricarboxylic acid
"""
"""
Classifies: tricarboxylic acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tricarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a tricarboxylic acid based on its SMILES string.
    A tricarboxylic acid is an oxoacid containing exactly three carboxy groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tricarboxylic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens to ensure all hydrogens are considered
    mol = Chem.AddHs(mol)

    # Define SMARTS pattern for carboxylic acid group (protonated or deprotonated)
    carboxylic_acid_pattern = Chem.MolFromSmarts('[C;X3](=O)[O;H1,-1]')

    # Find all carboxylic acid groups
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    num_carboxylic_acids = len(carboxylic_acid_matches)

    # Checking for exactly three carboxylic acid groups
    if num_carboxylic_acids != 3:
        return False, f"Contains {num_carboxylic_acids} carboxylic acid group(s), expected exactly 3"

    # Ensure that the carboxyl groups are not part of esters, amides, or other derivatized forms
    # We can check for adjacent atoms that might indicate derivatization
    ester_amide_pattern = Chem.MolFromSmarts('[C;X3](=O)[O,N][C]')
    if mol.HasSubstructMatch(ester_amide_pattern):
        return False, "Carboxylic acid groups are derivatized (e.g., esters or amides)"

    # Optionally, check that the molecule is an oxoacid (contains at least one oxo group)
    oxo_group_pattern = Chem.MolFromSmarts('[OX1]=[C]')
    if not mol.HasSubstructMatch(oxo_group_pattern):
        return False, "No oxo group found, not an oxoacid"

    return True, "Contains exactly three carboxylic acid groups and is a tricarboxylic acid"

__metadata__ = {   'chemical_class': {   'id': None,
                              'name': 'tricarboxylic acid',
                              'definition': 'An oxoacid containing three carboxy groups.'},
        'config': {   'llm_model_name': 'lbl/claude-sonnet',
                      'f1_threshold': 0.8,
                      'max_attempts': 5,
                      'max_positive_instances': None,
                      'max_positive_to_test': None,
                      'max_negative_to_test': None,
                      'max_positive_in_prompt': 50,
                      'max_negative_in_prompt': 20,
                      'max_instances_in_prompt': 100,
                      'test_proportion': 0.1}}