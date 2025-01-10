"""
Classifies: CHEBI:37142 organoiodine compound
"""
"""
Classifies: CHEBI:24863 organoiodine compound
Definition: A compound containing at least one carbon-iodine bond.
"""
from rdkit import Chem

def is_organoiodine_compound(smiles: str):
    """
    Determines if a molecule is an organoiodine compound based on its SMILES string.
    An organoiodine compound must contain at least one carbon-iodine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organoiodine compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for C-I bonds using SMARTS pattern
    # This will match any carbon-iodine bond regardless of charge/isotope
    ci_pattern = Chem.MolFromSmarts('C-I')
    
    if mol.HasSubstructMatch(ci_pattern):
        # Count how many C-I bonds for the explanation
        matches = len(mol.GetSubstructMatches(ci_pattern))
        return True, f"Contains {matches} carbon-iodine bond{'s' if matches > 1 else ''}"
    
    # Alternative pattern to catch aromatic carbon-iodine bonds
    ci_aromatic_pattern = Chem.MolFromSmarts('c-I')
    
    if mol.HasSubstructMatch(ci_aromatic_pattern):
        matches = len(mol.GetSubstructMatches(ci_aromatic_pattern))
        return True, f"Contains {matches} aromatic carbon-iodine bond{'s' if matches > 1 else ''}"

    return False, "No carbon-iodine bonds found"