"""
Classifies: CHEBI:25409 monoterpenoid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.
    A monoterpenoid typically has a carbon skeleton derived from a monoterpene
    with around 10 carbons, potentially rearranged or modified.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Check for approximately 10 carbon atoms
    if c_count < 8 or c_count > 12:
        return False, f"Carbon count {c_count} not typical for monoterpenoid"
    
    # We can't definitively determine rearrangements or skeletal removals
    # without specific patterns, so we'll return (None, None) to indicate uncertainty
    return None, "Unable to definitively classify due to complexity of rearrangement identification"