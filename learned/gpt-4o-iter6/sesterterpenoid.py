"""
Classifies: CHEBI:26660 sesterterpenoid
"""
from rdkit import Chem

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.

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

    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    # Adjust range to consider unusual configurations or modifications
    if c_count < 20 or c_count > 40:
        return False, f"Carbon count ({c_count}) not typical or too large for a sesterterpenoid"

    # Check for the presence of isoprene and related patterns
    # Consider rearranged/demethylated versions typical in sesterterpenoids
    isoprene_patterns = [
        Chem.MolFromSmarts("C=C(C)C"),  # Basic isoprene unit
        Chem.MolFromSmarts("C=CC(C)C"), # Rearranged form
        Chem.MolFromSmarts("CC(C)=C"),  # Different alkene location
    ]
    
    # Check if any isoprene or similar pattern is present
    if not any(mol.HasSubstructMatch(pattern) for pattern in isoprene_patterns):
        return False, "No isoprene-like units found, which are common in terpenoids"

    # Further terpenoid-specific structures could be checked here, e.g., certain cyclic structures

    # If structure fits broadened criteria
    return True, "Likely a sesterterpenoid based on carbon count and structural features"