"""
Classifies: CHEBI:4194 D-hexose
"""
from rdkit import Chem

def is_D_hexose(smiles: str):
    """
    Determines if a molecule is a D-hexose based on its SMILES string.
    A D-hexose is a six-carbon monosaccharide with a specific stereochemistry.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a D-hexose, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check for hexose structure (six carbons with appropriate hydroxyl groups)
    # D-hexoses have at least one aldehyde group or form a hemiacetal/pyranose
    hexose_pattern = Chem.MolFromSmarts("C(C(C(C(C(C=O|O1)))O)O)O")
    if not mol.HasSubstructMatch(hexose_pattern):
        return False, "Does not match hexose pattern"
    
    # Check for the presence of chirality at four chiral centers typical for D-hexoses
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if len(chiral_centers) < 4:
        return False, f"Expected at least four chiral centers, found {len(chiral_centers)}"
    
    # Further checks for stereochemistry (e.g., at specific carbons for exact D configuration)
    # This is a simplified version, a more accurate match would require specific D-hexose stereochemistry
    # Check for at least two `@` stereochemical markers
    if smiles.count('@') < 2:
        return False, "Not enough stereochemical markers for D-hexose"

    return True, "Contains D-hexose pattern with sufficient stereochemical detail"