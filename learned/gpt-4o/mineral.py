"""
Classifies: CHEBI:46662 mineral
"""
from rdkit import Chem

def is_mineral(smiles: str):
    """
    Determines if a molecule is a mineral based on its SMILES string.
    Minerals typically contain metal ions, sulfides, oxides, phosphates, or are metalloids.

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
    ]

    # Check for presence of metal ions or metalloids
    contains_metal_or_metalloid = any(atom.GetAtomicNum() in metal_and_metalloid_atomic_numbers for atom in mol.GetAtoms())
    if not contains_metal_or_metalloid:
        return False, "No metal ions or metalloids found in structure"
    
    # Check for presence of key anions and functional groups typical in minerals
    mineral_group_smarts = [
        '[O-]S(=O)(=O)[O-]',  # sulfate
        'C(=O)([O-])[O-]',    # carbonate
        'P(=O)([O-])([O-])[O-]',  # phosphate
        '[N+]([O-])=O'        # nitrate
    ]
    
    for group_smarts in mineral_group_smarts:
        group_pattern = Chem.MolFromSmarts(group_smarts)
        if mol.HasSubstructMatch(group_pattern):
            return True, "Contains metal ions/metalloids and common mineral anion"
    
    # Consider sulfides and oxides common in minerals
    sulfide_oxide_patterns = [
        '[S]',  # sulfur presence indicating possible sulfide
        '[O]'    # oxygen presence - possible oxides or hydrates
    ]
    
    for pattern in sulfide_oxide_patterns:
        substruct_pattern = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(substruct_pattern):
            return True, "Contains common mineral-forming elements like S or O"
    
    return False, "Does not match typical mineral structure"