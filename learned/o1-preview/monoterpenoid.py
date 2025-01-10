"""
Classifies: CHEBI:25409 monoterpenoid
"""
"""
Classifies: monoterpenoid
"""
from rdkit import Chem

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.
    A monoterpenoid is a terpenoid derived from a monoterpene (C10 skeleton).
    The term includes compounds in which the C10 skeleton of the parent monoterpene
    has been rearranged or modified by the removal of one or more skeletal atoms
    (generally methyl groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Monoterpenoids are derived from monoterpenes (C10 skeleton)
    # Allow for modifications, so require carbons between 7 and 12
    if c_count < 7 or c_count > 12:
        return False, f"Number of carbon atoms ({c_count}) not in range typical of monoterpenoids (7-12)"
    
    # Check for presence of oxygen atoms (typical for terpenoids)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 1:
        return False, "No oxygen atoms found; monoterpenoids typically contain oxygen functional groups"

    return True, "Molecule has characteristics of a monoterpenoid (appropriate carbon count and contains oxygen)"