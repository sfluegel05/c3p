"""
Classifies: CHEBI:72544 flavonoids
"""
"""
Classifies: CHEBI:47916 flavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_flavonoids(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    Flavonoids are characterized by a C6-C3-C6 skeleton, often with hydroxyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the C6-C3-C6 skeleton pattern
    # This pattern matches two aromatic rings connected by a three-carbon chain, allowing for variations
    skeleton_pattern = Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX4]-[CX4]-[CX4]-[c]2[c][c][c][c][c]2")
    skeleton_pattern_alt1 = Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX4]-[CX4]-[c]2[c][c][c][c][c]2")
    skeleton_pattern_alt2 = Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX4]-[c]2[c][c][c][c][c]2")
    skeleton_pattern_alt3 = Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX4]-[CX4]-[CX4]-[c]2[c][c][c][c][c]2")
    skeleton_pattern_alt4 = Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX4]-[CX4]-[c]2[c][c][c][c][c]2")
    skeleton_pattern_alt5 = Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX4]-[c]2[c][c][c][c][c]2")

    if not (mol.HasSubstructMatch(skeleton_pattern) or 
            mol.HasSubstructMatch(skeleton_pattern_alt1) or 
            mol.HasSubstructMatch(skeleton_pattern_alt2) or 
            mol.HasSubstructMatch(skeleton_pattern_alt3) or 
            mol.HasSubstructMatch(skeleton_pattern_alt4) or 
            mol.HasSubstructMatch(skeleton_pattern_alt5)):
        return False, "No C6-C3-C6 skeleton found"

    # Count aromatic rings to ensure there are at least two
    aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if aromatic_rings < 2:
        return False, f"Found {aromatic_rings} aromatic rings, need at least 2"

    # Check for hydroxyl groups attached to the aromatic rings
    hydroxyl_pattern = Chem.MolFromSmarts("[c][OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 1:
        return False, f"Found {len(hydroxyl_matches)} hydroxyl groups, need at least 1"

    # Check molecular weight - flavonoids typically >200 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for flavonoid"

    return True, "Contains C6-C3-C6 skeleton with at least one hydroxyl group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:47916',
                          'name': 'flavonoid',
                          'definition': 'Any organic molecular entity whose structure is based on derivatives of a phenyl-substituted 1-phenylpropane possessing a C15 or C16 skeleton, or such a structure which is condensed with a C6-C3 lignan precursors.',
                          'parents': ['CHEBI:24431', 'CHEBI:50753']},
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