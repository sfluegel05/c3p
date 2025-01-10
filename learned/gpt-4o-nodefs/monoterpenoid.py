"""
Classifies: CHEBI:25409 monoterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.
    A monoterpenoid is a terpene consisting of two isoprene units, typically containing 10 carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a monoterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 10:
        return False, f"Expected 10 carbons, found {c_count}"
    
    # Check presence of isoprene unit patterns where appropriate
    # Isoprene C5H8 pattern: C=C-C-C=C
    isoprene_pattern = Chem.MolFromSmarts("C=C-C-C=C")
    if not mol.HasSubstructMatch(isoprene_pattern):
        return False, "No isoprene-like patterns found"
    
    # Confirm isoprene-derived structure availability
    # Check for typical cyclic structures in monoterpenoids (e.g., cyclohexane)
    cyclic_structure = Chem.MolFromSmarts("C1CCCCC1")
    if mol.HasSubstructMatch(cyclic_structure):
        return True, "Monoterpenoid with cyclic structure recognized"
    
    # Additionally, check for common moieties in monoterpenoids
    common_motif_1 = Chem.MolFromSmarts("C1=CC=CC=C1")  # Example for aromatic cycles
    common_motif_2 = Chem.MolFromSmarts("C=C(C)C")  # Example of patterns from terpenes

    if mol.HasSubstructMatch(common_motif_1) or mol.HasSubstructMatch(common_motif_2):
        return True, "Monoterpenoid structure with common moiety recognized"
    
    # Fallback reason if matches are inconclusive
    return False, "Structure differs from typical monoterpenoids"