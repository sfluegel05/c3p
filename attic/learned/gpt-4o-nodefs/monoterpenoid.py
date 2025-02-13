"""
Classifies: CHEBI:25409 monoterpenoid
"""
from rdkit import Chem

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.
    Monoterpenoids typically consist of two isoprene units often totaling 10 carbons,
    but with structural diversities such as cyclic, acyclic and functionalized versions.
    
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

    # Count carbon atoms, increase range to account for diversity
    c_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    c_count = len(c_atoms)
    if c_count < 8 or c_count > 20:
        return False, f"Carbon count {c_count} is atypical for a monoterpenoid"

    # Recognize cyclic patterns
    cyclic_patterns = [
        Chem.MolFromSmarts("C1CCCCC1"),  # Simple cyclohexane
        Chem.MolFromSmarts("C1CC2CCC(C2)C1"),  # Bicyclic example
        Chem.MolFromSmarts("C1=CC2CCCCC2C=C1"),  # Substituted cyclic motifs
    ]
    for pattern in cyclic_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Monoterpenoid with recognized cyclic structure identified"

    # Expanded acyclic structure patterns
    acyclic_patterns = [Chem.MolFromSmarts(pattern) for pattern in [
        "C=C(C)C",       # Isoprene unit
        "CCC=C",         # Extended chain
        "CC(C)=C",       # Other acyclic derivative
        "C=C(CC)CC",     # Typical acyclic monoterpene
    ]]
    for pattern in acyclic_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Monoterpenoid with known acyclic structure identified"
    
    # Check for characteristic functional group substructures, including more diversity
    functional_groups = [Chem.MolFromSmarts(smarts) for smarts in [
        "O", "C(=O)O",   # Alcohols and esters
        "C=O",           # Ketones
        "C1=CC=CC=C1",   # Aromatic
    ]]
    for group in functional_groups:
        if mol.HasSubstructMatch(group):
            return True, "Monoterpenoid with typical functional group identified"

    # Extra consideration for aromaticity
    aromatic_ring = Chem.MolFromSmarts("c1ccccc1")
    if mol.HasSubstructMatch(aromatic_ring):
        return True, "Monoterpenoid with aromatic structure detected"

    # Default return if no clear monoterpenoid patterns or structures identified
    return False, "Does not match typical monoterpenoid patterns or structures"