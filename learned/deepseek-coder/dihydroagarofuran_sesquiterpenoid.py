"""
Classifies: CHEBI:71548 dihydroagarofuran sesquiterpenoid
"""
"""
Classifies: CHEBI:XXXXX dihydroagarofuran sesquiterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_dihydroagarofuran_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a dihydroagarofuran sesquiterpenoid based on its SMILES string.
    A dihydroagarofuran sesquiterpenoid is a sesquiterpenoid with a dihydroagarofuran skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dihydroagarofuran sesquiterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core dihydroagarofuran skeleton pattern
    dihydroagarofuran_pattern = Chem.MolFromSmarts("[C@H]12[C@@H](OC(=O)*)[C@H](OC(=O)*)[C@]3(COC(=O)*)[C@H](OC(=O)*)[C@@H](OC(=O)*)[C@]4([H])[C@@H](OC(=O)*)[C@]3(O[C@@]14C)[C@@]2(C)O")
    if not mol.HasSubstructMatch(dihydroagarofuran_pattern):
        return False, "No dihydroagarofuran skeleton found"

    # Check for sesquiterpenoid characteristics (15 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 15:
        return False, f"Expected 15 carbons for sesquiterpenoid, found {c_count}"

    # Check for oxygen-containing functional groups (esters, hydroxyls)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 3:
        return False, f"Expected at least 3 ester groups, found {len(ester_matches)}"

    # Check for the presence of a fused tricyclic system
    ring_info = mol.GetRingInfo()
    if len(ring_info.AtomRings()) < 3:
        return False, "Expected at least 3 rings in the fused tricyclic system"

    return True, "Contains dihydroagarofuran skeleton with sesquiterpenoid characteristics"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:XXXXX',
                          'name': 'dihydroagarofuran sesquiterpenoid',
                          'definition': 'Any sesquiterpenoid with a dihydroagarofuran skeleton.',
                          'parents': ['CHEBI:XXXXX', 'CHEBI:XXXXX']},
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