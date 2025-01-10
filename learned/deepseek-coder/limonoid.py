"""
Classifies: CHEBI:39434 limonoid
"""
"""
Classifies: CHEBI:26154 limonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_limonoid(smiles: str):
    """
    Determines if a molecule is a limonoid based on its SMILES string.
    Limonoids are highly oxygenated triterpenoids with a 4,4,8-trimethyl-17-furanylsteroid skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a limonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a furan ring (C=COC)
    furan_pattern = Chem.MolFromSmarts("c1ccoc1")
    if not mol.HasSubstructMatch(furan_pattern):
        return False, "No furan ring found"

    # Check for multiple oxygen-containing functional groups (hydroxyl, carbonyl, ester)
    oxygen_pattern = Chem.MolFromSmarts("[OX2H1,O=,O-C=O]")
    oxygen_matches = mol.GetSubstructMatches(oxygen_pattern)
    if len(oxygen_matches) < 3:
        return False, f"Found {len(oxygen_matches)} oxygen-containing groups, need at least 3"

    # Check for the triterpenoid skeleton (C30H48 or similar)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20 or c_count > 40:
        return False, f"Carbon count {c_count} is not typical for a triterpenoid"

    # Check for the presence of a steroid-like skeleton (4,4,8-trimethyl pattern)
    steroid_pattern = Chem.MolFromSmarts("[CX4H3][CX4H3][CX4H3]")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No 4,4,8-trimethyl pattern found"

    # Check molecular weight (limonoids typically >400 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for a limonoid"

    # Count rotatable bonds (limonoids typically have several)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Too few rotatable bonds for a limonoid"

    # Check for the presence of multiple rings (limonoids typically have several rings)
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    if n_rings < 4:
        return False, "Too few rings for a limonoid"

    return True, "Contains furan ring, multiple oxygen-containing groups, and a triterpenoid skeleton"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26154',
                          'name': 'limonoid',
                          'definition': 'Any triterpenoid that is highly oxygenated and has a prototypical structure either containing or derived from a precursor with a 4,4,8-trimethyl-17-furanylsteroid skeleton.',
                          'parents': ['CHEBI:26153', 'CHEBI:26152']},
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