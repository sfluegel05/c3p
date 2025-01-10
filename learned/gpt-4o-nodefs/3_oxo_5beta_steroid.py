"""
Classifies: CHEBI:1624 3-oxo-5beta-steroid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_5beta_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5beta-steroid based on its SMILES string.
    A 3-oxo-5beta-steroid should have a steroid backbone with a ketone group at the 3-position
    and specific stereochemistry denoted as '5beta'.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 3-oxo-5beta-steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 3-oxo group pattern (ketone at 3-position)
    oxo_pattern_3 = Chem.MolFromSmarts("C(=O)[C@@H]")  # Assuming this pattern for '3-oxo'
    if not mol.HasSubstructMatch(oxo_pattern_3):
        return False, "3-oxo group not found"
    
    # Stereochemistry at 5beta: This would be part of a complete pattern for a steroid with 5beta annotation
    # Complex, assumed 'beta' refers configuration that requires more detailed chiral recognition
    # Simple check for presence of chiral centers
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    
    # Look for the steroid scaffold pattern (ABCD ring system)
    steroid_pattern = Chem.MolFromSmarts("C1C2C3C4") # Extremely simplified steroid nucleus
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Steroid backbone not found"
    
    # Optional: More detailed chiral pattern matching
    has_beta_specificity = any('beta' in center for center in chiral_centers)
    if not has_beta_specificity:
        return False, f"5beta stereochemistry not found; chiral centers: {chiral_centers}"
    
    # Additional pattern specifics for merged features could be added here as needed

    return True, "3-oxo-5beta-steroid characteristic patterns found"

# The patterns in the code are placeholders. More precise pattern definitions
# are required to handle actual stereochemical configurations.