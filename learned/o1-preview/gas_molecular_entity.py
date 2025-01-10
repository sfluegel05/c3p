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

    This function checks if the molecule consists only of main group elements and estimates its physical state at STP based on molecular weight.

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
    main_group_atomic_numbers = {
        # Elements of groups 1, 2, and 13-18
        1, 2,                                    # H, He
        3, 4,                                    # Li, Be
        5, 6, 7, 8, 9, 10,                       # B, C, N, O, F, Ne
        11, 12,                                  # Na, Mg
        13, 14, 15, 16, 17, 18,                  # Al, Si, P, S, Cl, Ar
        31, 32, 33, 34, 35, 36,                  # Ga, Ge, As, Se, Br, Kr
        49, 50, 51, 52, 53, 54,                  # In, Sn, Sb, Te, I, Xe
        81, 82, 83, 84, 85, 86                   # Tl, Pb, Bi, Po, At, Rn
    }

    # Check that all atoms are main group elements
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in main_group_atomic_numbers:
            elem_symbol = atom.GetSymbol()
            return False, f"Contains non-main group element: {elem_symbol} (atomic number {atomic_num})"

    # Calculate molecular weight
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)

    # Define molecular weight threshold (approximate)
    weight_threshold = 60.0  # g/mol

    # Special handling for noble gases (monoatomic molecules)
    noble_gases = {2, 10, 18, 36, 54, 86}  # He, Ne, Ar, Kr, Xe, Rn
    if mol.GetNumAtoms() == 1:
        atomic_num = mol.GetAtomWithIdx(0).GetAtomicNum()
        if atomic_num in noble_gases:
            return True, f"Single atom of noble gas element: {mol.GetAtomWithIdx(0).GetSymbol()}"

    # Classify based on molecular weight
    if mol_weight <= weight_threshold:
        return True, f"Molecular weight {mol_weight:.2f} g/mol is below threshold of {weight_threshold} g/mol, likely a gas at STP"
    else:
        return False, f"Molecular weight {mol_weight:.2f} g/mol exceeds threshold of {weight_threshold} g/mol, likely not a gas at STP"

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
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None
}