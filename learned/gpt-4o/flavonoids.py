"""
Classifies: CHEBI:72544 flavonoids
"""
from rdkit import Chem

def is_flavonoids(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for a generic flavonoid core
    flavonoid_core_pattern = Chem.MolFromSmarts("c1ccc2c(c1)C(=O)c3c(cccc3)c2")  # General flavone core
    if flavonoid_core_pattern is None:
        return False, "Critical error: Could not compile flavonoid core pattern"

    if not mol.HasSubstructMatch(flavonoid_core_pattern):
        return False, "No flavonoid core found"

    # Check for at least two aromatic rings
    aromatic_rings = [ring for ring in mol.GetRingInfo().AtomRings()
                      if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)]
    if len(aromatic_rings) < 2:
        return False, "Too few aromatic rings for flavonoid classification"

    # Check for common functional groups (e.g., hydroxyl)
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    if hydroxyl_pattern is None:
        return False, "Critical error: Could not compile hydroxyl pattern"

    hydroxyl_count = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if hydroxyl_count == 0:
        return False, "No hydroxyl group found which is typical for flavonoids"

    return True, "Flavonoid detected with a recognized core and functional group presence"