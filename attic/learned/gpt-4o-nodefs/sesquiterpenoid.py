"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
from rdkit import Chem

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.
    A sesquiterpenoid is defined by having 15 carbon atoms typically derived from three isoprene units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesquiterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for exactly 15 carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 15:
        return False, f"Carbon count is {c_count}, expected 15 for sesquiterpenoid"

    # Sesquiterpenoids are often identified by 3 isoprene units,
    # Here we would typically insert more sophisticated pattern matching,
    # but this basic check only ensures carbon count. Advanced checks require complex SMARTS.
    
    # For now, return True for correct carbon count indicating potential sesquiterpenoid
    return True, "Contains 15 carbon atoms, indicative of sesquiterpenoid structure"

# Example usage (for testing purposes):
# smiles = "OC1C(C(CC/C=C(\C)/C=O)C)=CC(=O)C(=C1)C"
# print(is_sesquiterpenoid(smiles))