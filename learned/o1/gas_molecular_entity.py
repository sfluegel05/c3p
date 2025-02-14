"""
Classifies: CHEBI:138675 gas molecular entity
"""
# gas_molecular_entity.py
"""
Classifies: CHEBI gas molecular entity
"""
from rdkit import Chem

def is_gas_molecular_entity(smiles: str):
    """
    Determines if a molecule is a main group molecular entity that is gaseous at STP.
    Uses Joback's method to estimate the boiling point of the molecule.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a gas molecular entity, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for net charge; exclude ions
    if Chem.GetFormalCharge(mol) != 0:
        return False, "Molecule is an ion with non-zero net charge, unlikely to be gaseous at STP"
    
    # Define Joback group contributions
    joback_groups = [
        {'name': 'CH3', 'smarts': '[CX4H3]', 'delta_tb': 23.58},
        {'name': 'CH2', 'smarts': '[CX4H2]', 'delta_tb': 22.88},
        {'name': 'CH',  'smarts': '[CX4H1]', 'delta_tb': 21.74},
        {'name': 'C',   'smarts': '[CX4H0]', 'delta_tb': 20.15},
        {'name': 'CH3 (aromatic)', 'smarts': '[cH3]', 'delta_tb': 17.21},
        {'name': 'CH2 (aromatic)', 'smarts': '[cH2]', 'delta_tb': 16.33},
        {'name': 'CH (aromatic)',  'smarts': '[cH1]', 'delta_tb': 13.60},
        {'name': 'C (aromatic)',   'smarts': '[cH0]', 'delta_tb': 12.65},
        {'name': 'C=C', 'smarts': 'C=C', 'delta_tb': 1.10},
        {'name': 'C#C', 'smarts': 'C#C', 'delta_tb': -27.87},
        {'name': 'C=O (aldehyde)', 'smarts': '[CX3H1](=O)[#6]', 'delta_tb': 26.00},
        {'name': 'C=O (ketone)', 'smarts': '[#6][CX3](=O)[#6]', 'delta_tb': 20.00},
        {'name': 'OH (alcohol)', 'smarts': '[OX2H]', 'delta_tb': 42.57},
        {'name': 'O (ether)', 'smarts': '[OD2]([#6])[#6]', 'delta_tb': 17.85},
        {'name': 'N (amine)', 'smarts': '[NX3H2]', 'delta_tb': 34.34},
        {'name': 'N (secondary amine)', 'smarts': '[NX3H1]([#6])[#6]', 'delta_tb': 30.50},
        {'name': 'N (tertiary amine)', 'smarts': '[NX3]([#6])([#6])[#6]', 'delta_tb': 27.00},
        {'name': 'F', 'smarts': '[F]', 'delta_tb': -15.92},
        {'name': 'Cl', 'smarts': '[Cl]', 'delta_tb': 9.17},
        {'name': 'Br', 'smarts': '[Br]', 'delta_tb': 22.80},
        {'name': 'I', 'smarts': '[I]', 'delta_tb': 36.37},
        # Add more functional groups as needed
    ]
    
    # Initialize total contribution
    delta_tb_total = 0.0
    
    # Create a copy of the molecule for group assignment
    mol_copy = Chem.AddHs(mol)
    
    # Create an atom mapping to avoid double counting
    atom_marks = [0] * mol_copy.GetNumAtoms()
    
    # Iterate over Joback groups and count occurrences
    for group in joback_groups:
        pattern = Chem.MolFromSmarts(group['smarts'])
        matches = mol_copy.GetSubstructMatches(pattern)
        for match in matches:
            # Check if any atom in the match is already marked
            if any(atom_marks[atom_idx] for atom_idx in match):
                continue  # Skip if already counted
            # Mark atoms as counted
            for atom_idx in match:
                atom_marks[atom_idx] = 1
            # Add group contribution
            delta_tb_total += group['delta_tb']
    
    # Estimate the boiling point
    estimated_tb = 198.0 + delta_tb_total  # Boiling point in Celsius
    
    # Determine if the molecule is gaseous at STP
    if estimated_tb <= 0.0:
        return True, f"Estimated boiling point {estimated_tb:.2f}°C suggests the molecule is gaseous at STP"
    else:
        return False, f"Estimated boiling point {estimated_tb:.2f}°C suggests the molecule is not gaseous at STP"