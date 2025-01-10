"""
Classifies: CHEBI:25409 monoterpenoid
"""
from rdkit import Chem

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.
    Monoterpenoids typically consist of 10 carbon atoms coming from two isoprene units.
    They might be cyclic or acyclic, and exhibit functional groups.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a monoterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES with RDKit
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbon atoms, allow slight variation
    c_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    c_count = len(c_atoms)
    if c_count < 8 or c_count > 16:
        return False, f"Carbon count {c_count} is atypical for a monoterpenoid"

    # Recognize cyclic patterns often found in monoterpenoids
    cyclic_patterns = [Chem.MolFromSmarts("C1CCCCC1"), Chem.MolFromSmarts("C1CC2CCC(C2)C1")]
    for pattern in cyclic_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Monoterpenoid with recognized cyclic structure identified"

    # Acyclic structure patterns that diverge slightly from a purely linear isoprene base
    acyclic_patterns = [Chem.MolFromSmarts(pattern) for pattern in [
        "C=C(C)C",  # Single isoprene unit
        "CCC=C",    # Extended chain structures
        "CC(C)=C"   # Other derivative motifs
    ]]
    
    for pattern in acyclic_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Monoterpenoid with known acyclic structure identified"
    
    # Check for characteristic functional group substructures - e.g., alcohols, esters
    functional_groups = [Chem.MolFromSmarts(smarts) for smarts in ["O", "C(=O)O"]]
    for group in functional_groups:
        if mol.HasSubstructMatch(group):
            return True, "Monoterpenoid with typical functional group identified"

    # Check for aromatic characteristics
    aromatic_ring = Chem.MolFromSmarts("c1ccccc1")
    if mol.HasSubstructMatch(aromatic_ring):
        return True, "Monoterpenoid with aromatic structure detected"

    # Default return if no clear monoterpenoid patterns or structures were identified
    return False, "Does not match typical monoterpenoid patterns or structures"