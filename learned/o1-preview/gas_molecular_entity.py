"""
Classifies: CHEBI:138675 gas molecular entity
"""
"""
Classifies: gas molecular entity
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_gas_molecular_entity(smiles: str):
    """
    Determines if a molecule is a gas molecular entity based on its SMILES string.
    A gas molecular entity is any main group molecular entity that is gaseous at standard temperature and pressure (STP; 0°C and 100 kPa).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a gas molecular entity, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define main group elements (atomic numbers)
    main_group_atomic_numbers = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                                 11, 12, 13, 14, 15, 16, 17, 18,
                                 19, 20, 31, 32, 33, 34, 35, 36,
                                 37, 38, 49, 50, 51, 52, 53, 54,
                                 55, 56, 81, 82, 83, 84, 85, 86,
                                 87, 88}

    # Check that all atoms are main group elements
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in main_group_atomic_numbers:
            elem_symbol = atom.GetSymbol()
            return False, f"Contains non-main group element: {elem_symbol} (atomic number {atomic_num})"

    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)

    # Set molecular weight threshold (e.g., 85 g/mol)
    mol_wt_threshold = 85.0
    if mol_wt > mol_wt_threshold:
        return False, f"Molecular weight {mol_wt:.2f} g/mol exceeds threshold of {mol_wt_threshold} g/mol"

    # Assume that molecule is gas if below threshold
    return True, f"Molecule is likely a gas at STP (molecular weight {mol_wt:.2f} g/mol)"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:22318',
        'name': 'gas molecular entity',
        'definition': 'Any main group molecular entity that is gaseous at standard temperature and pressure (STP; 0°C and 100 kPa).',
        'parents': ['CHEBI:33565']
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
    'stdout': None
}