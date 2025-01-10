"""
Classifies: CHEBI:46662 mineral
"""
"""
Classifies: CHEBI:46662 mineral
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_mineral(smiles: str):
    """
    Determines if a molecule is a mineral based on its SMILES string.
    A mineral is generally a chemical substance that is normally crystalline and formed through geological processes.
    This includes inorganic compounds, salts, and certain naturally occurring amorphous substances.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mineral, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of metal ions (common in minerals)
    metal_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() in [3, 11, 12, 13, 19, 20, 22, 24, 25, 26, 27, 28, 29, 30, 38, 40, 42, 47, 50, 56, 74, 78, 79, 80, 82, 83, 88, 92]]
    if not metal_atoms:
        return False, "No metal ions found, which are common in minerals"

    # Check for the presence of anions (common in minerals)
    anions = [atom for atom in mol.GetAtoms() if atom.GetFormalCharge() < 0]
    if not anions:
        return False, "No anions found, which are common in minerals"

    # Check for common mineral patterns (e.g., sulfates, carbonates, phosphates)
    sulfate_pattern = Chem.MolFromSmarts("[O-]S(=O)(=O)[O-]")
    carbonate_pattern = Chem.MolFromSmarts("[O-]C(=O)[O-]")
    phosphate_pattern = Chem.MolFromSmarts("[O-]P(=O)([O-])[O-]")
    silicate_pattern = Chem.MolFromSmarts("[O-][Si]([O-])([O-])[O-]")
    
    if mol.HasSubstructMatch(sulfate_pattern):
        return True, "Contains sulfate group, common in minerals"
    if mol.HasSubstructMatch(carbonate_pattern):
        return True, "Contains carbonate group, common in minerals"
    if mol.HasSubstructMatch(phosphate_pattern):
        return True, "Contains phosphate group, common in minerals"
    if mol.HasSubstructMatch(silicate_pattern):
        return True, "Contains silicate group, common in minerals"

    # Check for simple salts (e.g., NaCl, KCl)
    if len(mol.GetAtoms()) <= 3:
        return True, "Simple salt, common in minerals"

    # Check for hydrates (common in minerals)
    hydrate_pattern = Chem.MolFromSmarts("[OH2]")
    if mol.HasSubstructMatch(hydrate_pattern):
        return True, "Contains water of hydration, common in minerals"

    # If none of the above patterns match, return False
    return False, "Does not match common mineral patterns"

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:46662',
                          'name': 'mineral',
                          'definition': 'In general, a mineral is a chemical substance that is normally crystalline formed and has been formed as a result of geological processes. The term also includes metamict substances (naturally occurring, formerly crystalline substances whose crystallinity has been destroyed by ionising radiation) and can include naturally occurring amorphous substances that have never been crystalline ('mineraloids') such as georgite and calciouranoite as well as substances formed by the action of geological processes on bigenic compounds ('biogenic minerals').'},
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