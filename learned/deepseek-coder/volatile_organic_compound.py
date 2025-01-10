"""
Classifies: CHEBI:134179 volatile organic compound
"""
"""
Classifies: CHEBI:17741 volatile organic compound
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_volatile_organic_compound(smiles: str):
    """
    Determines if a molecule is a volatile organic compound (VOC) based on its SMILES string.
    A VOC is defined as any organic compound with an initial boiling point ≤ 250°C at 101.3 kPa.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a VOC, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule is organic (contains carbon)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count == 0:
        return False, "Not an organic compound (no carbon atoms)"

    # Calculate molecular weight
    mol_wt = Descriptors.MolWt(mol)

    # Check for molecular weight as a proxy for boiling point
    # VOCs typically have lower molecular weights
    if mol_wt > 250:  # Empirical threshold based on typical VOC molecular weights
        return False, f"Molecular weight ({mol_wt:.2f}) too high for VOC"

    # Check for heavy atoms (non-hydrogen atoms)
    heavy_atom_count = mol.GetNumHeavyAtoms()
    if heavy_atom_count > 20:  # Empirical threshold for VOCs
        return False, f"Too many heavy atoms ({heavy_atom_count}) for VOC"

    # Check for functional groups that might increase boiling point
    # (e.g., multiple hydroxyl groups, long carbon chains)
    hydroxyl_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[OX2H]")))
    if hydroxyl_count > 1:
        return False, f"Too many hydroxyl groups ({hydroxyl_count}) for VOC"

    # Check for long carbon chains
    long_chain_pattern = Chem.MolFromSmarts("[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]")
    if mol.HasSubstructMatch(long_chain_pattern):
        return False, "Long carbon chain detected, likely not a VOC"

    # If all checks pass, classify as VOC
    return True, f"Molecular weight ({mol_wt:.2f}) and structure consistent with VOC criteria"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:17741',
                          'name': 'volatile organic compound',
                          'definition': 'Any organic compound having an initial '
                                        'boiling point less than or equal to '
                                        '250 degreeC (482 degreeF) measured at '
                                        'a standard atmospheric pressure of '
                                        '101.3 kPa.',
                          'parents': ['CHEBI:76971', 'CHEBI:76972']},
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