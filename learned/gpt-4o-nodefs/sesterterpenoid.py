"""
Classifies: CHEBI:26660 sesterterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.
    Sesterterpenoids are terpenoids consisting of 25 carbon atoms (five isoprene units),
    often with complex ring structures and various functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesterterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbon atoms to check if it has around 25 carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (20 <= c_count <= 30):
        return False, "Does not have 25 carbon atoms typical for sesterterpenoids"
    
    # Search for specific functional groups or ring structures, if possible
    # (This is a placeholder step. Specific substructures would need additional SMILES patterns.)
    # Example: Check for at least one ring structure
    if not mol.GetRingInfo().NumRings() > 0:
        return False, "No ring structure found, typically sesterterpenoids have complex rings"
    
    # Check for typical terpenoid-like features (isoprene units, common functional groups)
    # Note: This requires complex SMARTS patterns which are skipped for simplicity.
    
    # Generally, if the carbon count is close and has ring structures, posit it as a sesterterpenoid
    return True, "Carbon count and ring presence suggest a sesterterpenoid"

# Example usage:
# print(is_sesterterpenoid("O=C1C=C(C)[C@@H]...")) # Replace `...` with a complete SMILES from the example list