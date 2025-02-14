"""
Classifies: CHEBI:46662 mineral
"""
from rdkit import Chem

def is_mineral(smiles: str):
    """
    Determines if a molecule is a mineral based on its SMILES string.
    Minerals typically contain metal ions/metalloids and anions such as sulfates, phosphates, chlorides, oxides, or carbonates.

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
    
    # List of common metal ions and metalloids (atomic numbers for metals and some metalloids)
    metal_and_metalloid_atomic_numbers = [
        3, 11, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 37, 38, 
        39, 40, 42, 47, 48, 56, 57, 58, 79, 80, 81, 82, 83, 51, 33 # added Sb and As
        # Additional elements commonly found in minerals
    ]

    # Check for presence of metal ions or metalloids
    contains_metal_or_metalloid = any(atom.GetAtomicNum() in metal_and_metalloid_atomic_numbers for atom in mol.GetAtoms())
    if not contains_metal_or_metalloid:
        return False, "No metal ions or metalloids found in structure"
    
    # Check for presence of applicable anions in minerals
    mineral_group_smarts = [
        '[O-]S(=O)(=O)[O-]',  # sulfate
        'C(=O)([O-])[O-]',    # carbonate
        'P(=O)([O-])([O-])[O-]',  # phosphate
        '[N+]([O-])=O',       # nitrate
        'Cl',  # chloride common in many minerals
        '[O]',  # potential oxides/hydrates indicator
        '[S-]'  # potential sulfide indicator
    ]
    
    for group_smarts in mineral_group_smarts:
        group_pattern = Chem.MolFromSmarts(group_smarts)
        if mol.HasSubstructMatch(group_pattern):
            if '[C]' not in smiles or len(smiles) < 45: # Ensure it's not a large organic compound
                return True, "Contains metal ions/metalloids and suitable anion."

    return False, "Does not match typical mineral structure"