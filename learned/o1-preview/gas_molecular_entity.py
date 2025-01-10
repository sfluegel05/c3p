"""
Classifies: CHEBI:138675 gas molecular entity
"""
"""
Classifies: gas molecular entity
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors

def is_gas_molecular_entity(smiles: str):
    """
    Determines if a molecule is a gas molecular entity based on its SMILES string.
    A gas molecular entity is any main group molecular entity that is gaseous at standard temperature and pressure (STP; 0°C and 100 kPa).

    This function estimates the boiling point of the molecule using the Joback method.
    If the estimated boiling point is below 0°C (273.15 K), it is classified as a gas at STP.

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

    # Estimate boiling point using Joback method
    group_contributions = {
        # Group: Contribution to Tb (K)
        'CH3':   23.58,
        'CH2':   22.88,
        'CH':    21.74,
        'C':     20.73,
        'CH4':   24.44,
        'C=C':   1.00,
        'C#C':   -1.00,
        'C-O':   21.00,
        'O':     9.00,
        'OH':    26.00,
        'C=O':   13.00,
        'COOH':  48.50,
        'NH2':   16.00,
        'N':     14.00,
        'F':     9.00,
        'Cl':    15.00,
        'Br':    21.00,
        'I':     25.00,
        'S':     19.00,
        # Add more groups as needed
    }

    # Initialize sum of contributions
    Tb = 198.0  # Base temperature in Kelvin

    # Count functional groups
    functional_groups = {
        'CH3': Chem.MolFromSmarts('[CH3]'),
        'CH2': Chem.MolFromSmarts('[CH2]'),
        'CH':  Chem.MolFromSmarts('[CH]'),
        'C':   Chem.MolFromSmarts('[C]'),
        'C=C': Chem.MolFromSmarts('C=C'),
        'C#C': Chem.MolFromSmarts('C#C'),
        'OH':  Chem.MolFromSmarts('[OH]'),
        'NH2': Chem.MolFromSmarts('N([H])[H]'),
        'Cl':  Chem.MolFromSmarts('[Cl]'),
        'Br':  Chem.MolFromSmarts('[Br]'),
        'F':   Chem.MolFromSmarts('[F]'),
        'I':   Chem.MolFromSmarts('[I]'),
        # Add more functional groups as needed
    }

    for group_name, smarts in functional_groups.items():
        matches = mol.GetSubstructMatches(smarts)
        count = len(matches)
        contribution = group_contributions.get(group_name, 0.0)
        Tb += count * contribution
        # Debug print
        # print(f"Group {group_name}: count={count}, contribution={contribution}, total Tb={Tb}")

    # Special handling for noble gases and small molecules
    num_atoms = mol.GetNumAtoms()
    if num_atoms == 1:
        atomic_num = mol.GetAtomWithIdx(0).GetAtomicNum()
        noble_gases = {2, 10, 18, 36, 54, 86}  # He, Ne, Ar, Kr, Xe, Rn
        if atomic_num in noble_gases:
            Tb = Descriptors.AtomicNumber(mol) * 10  # Approximate boiling point
        else:
            Tb = 20.0  # Default for other single atoms

    # Classify based on estimated boiling point
    if Tb < 273.15:
        return True, f"Estimated boiling point {Tb:.2f} K is below 273.15 K (0°C), likely a gas at STP"
    else:
        return False, f"Estimated boiling point {Tb:.2f} K is above 273.15 K (0°C), not a gas at STP"

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
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None
}