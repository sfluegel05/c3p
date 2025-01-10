"""
Classifies: CHEBI:138675 gas molecular entity
"""
"""
Classifies: gas molecular entity
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_gas_molecular_entity(smiles: str):
    """
    Determines if a molecule is a gas molecular entity based on its SMILES string.
    A gas molecular entity is any main group molecular entity that is gaseous at standard temperature and pressure (STP; 0°C and 100 kPa).
    
    This function estimates the normal boiling point of the molecule using the Joback method.
    If the estimated boiling point is below 0°C (273.15 K), it is classified as a gas at STP.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is a gas molecular entity, False otherwise
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
    
    # Estimate boiling point using the Joback method
    try:
        tb = estimate_boiling_point_joback(mol)
    except Exception as e:
        return False, f"Error in boiling point estimation: {e}"
    
    # Classify based on boiling point
    if tb < 273.15:  # 0°C in Kelvin
        return True, f"Estimated boiling point {tb:.2f} K is below 273.15 K (0°C), likely a gas at STP"
    else:
        return False, f"Estimated boiling point {tb:.2f} K is above 273.15 K (0°C), likely not a gas at STP"
        
def estimate_boiling_point_joback(mol):
    """
    Estimates the normal boiling point (in Kelvin) of a molecule using the Joback method.
    
    Args:
        mol (Chem.Mol): RDKit molecule object
    
    Returns:
        float: Estimated boiling point in Kelvin
    """
    # Joback group contributions for boiling point
    group_contributions = {
        # Format: ('Smarts pattern', Tb increment in Kelvin)
        '[C;H3]':       23.58,  # -CH3
        '[C;H2]':       22.88,  # -CH2-
        '[C;H1]':       21.74,  # >CH-
        '[C;H0]':       20.56,  # >C<
        '[C;c]':        0.0,    # Aromatic carbons
        '[O;H1]':       9.09,   # -OH (alcohol)
        '[O;H0]':       9.04,   # >O (ether)
        '[N;H2]':       13.89,  # -NH2
        '[N;H1]':       13.60,  # >NH
        '[N;H0]':       14.13,  # >N-
        '[F]':          10.73,  # -F
        '[Cl]':         21.23,  # -Cl
        '[Br]':         21.86,  # -Br
        '[I]':          22.88,  # -I
        '[S;H0]':       19.22,  # >S<
        '[S;H1]':       13.81,  # -SH
        '[P]':          24.42,  # Phosphorus
        '[=O]':         7.97,   # C=O (carbonyl)
        '[#1]':         0.0,    # Hydrogen atoms (not counted)
    }
    
    # Initialize sum of group increments
    tb_increment = 198.0  # Base value in Kelvin
    
    # Create a copy of the molecule and add Hs to ensure all atoms are explicit
    mol = Chem.AddHs(mol)
    
    # Assign atom types based on hybridization and aromaticity
    for pattern, increment in group_contributions.items():
        smarts = Chem.MolFromSmarts(pattern)
        matches = mol.GetSubstructMatches(smarts)
        n_matches = len(matches)
        tb_increment += increment * n_matches

    return tb_increment

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
    'attempt': 3,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None
}