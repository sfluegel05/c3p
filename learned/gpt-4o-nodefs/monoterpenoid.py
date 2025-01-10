"""
Classifies: CHEBI:25409 monoterpenoid
"""
from rdkit import Chem

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.
    Monoterpenoids are terpenes consisting typically of 10 carbon atoms, but might vary slightly.
    They often contain two isoprene units, but structures can vary.
    
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
    
    # Count carbons, allow slight variation
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 9 or c_count > 15:
        return False, f"Carbon count {c_count} is atypical for a monoterpenoid"
    
    # Detect presence of cyclic and acyclic structures separately
    # Cyclic Example: detect typical cyclohexane or variant structures
    cyclic_pattern = Chem.MolFromSmarts("C1CCCCC1")
    if mol.HasSubstructMatch(cyclic_pattern):
        return True, "Monoterpenoid with cyclic structure identified"

    # Acyclic Patterns: Increased flexibility beyond simple isoprene motifs
    common_acyclic_motif = [Chem.MolFromSmarts(pattern) for pattern in [
        "C=C(C)C",  # Simple isoprenoid units
        "CC=C"      # Possible structure to detect with functional groups
    ]]

    for motif in common_acyclic_motif:
        if mol.HasSubstructMatch(motif):
            return True, "Monoterpenoid with known acyclic structure identified"
    
    # Check for presence of common aromatic or functional groups
    benzene_ring = Chem.MolFromSmarts("c1ccccc1")
    if mol.HasSubstructMatch(benzene_ring):
        return True, "Monoterpenoid with aromatic structure identified"

    return False, "Does not match typical monoterpenoid patterns or structures"