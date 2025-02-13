"""
Classifies: CHEBI:46662 mineral
"""
from rdkit import Chem

def is_mineral(smiles: str):
    """
    Determines if a molecule is a mineral based on its SMILES string.
    Minerals typically contain metal ions and can form crystalline structures.

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
    
    # List of common metal ions (atomic numbers for metals)
    metal_atomic_numbers = [3, 11, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 37, 38, 39, 40, 42, 47, 48, 50, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83]

    # Check for presence of metal ions
    contains_metal = any(atom.GetAtomicNum() in metal_atomic_numbers for atom in mol.GetAtoms())
    if not contains_metal:
        return False, "No metal ions found in structure"
    
    # Check for inorganic anions or functional groups
    inorganic_anions_smarts = [
        '[O-]S(=O)(=O)[O-]',  # sulfate
        'C(=O)([O-])[O-]',    # carbonate
        'P(=O)([O-])([O-])[O-]',  # phosphate
        '[N+]([O-])=O'        # nitrate
    ]
    
    for anion_smarts in inorganic_anions_smarts:
        anion_pattern = Chem.MolFromSmarts(anion_smarts)
        if mol.HasSubstructMatch(anion_pattern):
            return True, "Contains metal ions and common inorganic anion"
    
    # If no specific inorganic groups are found, but metal ion is present
    if contains_metal:
        return True, "Contains metal ions"

    return False, "Does not match typical mineral structure"