"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
"""
Classifies: CHEBI:24532 iridoid monoterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines if a molecule is an iridoid monoterpenoid based on its SMILES string.
    Iridoids typically have a cyclopentane ring fused to a six-membered oxygen heterocycle,
    or their seco-derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an iridoid monoterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic molecular checks
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150 or mol_wt > 800:
        return False, "Molecular weight outside typical range for iridoids"

    # Count atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 8 or c_count > 30:
        return False, "Carbon count outside typical range for iridoids"
    if o_count < 2:
        return False, "Insufficient oxygen atoms for iridoid structure"

    # Specific iridoid core patterns with more precise atom specifications
    core_patterns = [
        # Classic cyclopentane-pyran fusion with specific connectivity
        "[C;R2]1[C;R2][C;R2]2[C;R2]1[O;R1][C;R1][C;R1][C;R1][C;R2]2",
        # Cyclopentane-lactone fusion
        "[C;R2]1[C;R2][C;R2]2[C;R2]1[O;R1][C;R1](=O)[C;R1][C;R2]2",
        # Hemiacetal pattern with specific oxygen placement
        "[C;R2]1[C;R2][C;R2]2[C;R2]1[O;R1][C;R1]([O;H1])[C;R1][C;R2]2",
        # Common secoiridoid pattern with opened cyclopentane
        "[C;R0][C;R0](=O)[C;R0][C;R0][O;R1][C;R1]=[C;R1][C;R1]",
        # Specific pattern for genipin-type structures
        "[C;R2]1[C;R2][C;R2]2[O;R1][C;R1]=[C;R1][C;R2]2[C;R2]1",
    ]

    # Check for core skeleton
    is_core_found = False
    matched_pattern = None
    for pattern in core_patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat and mol.HasSubstructMatch(pat):
            is_core_found = True
            matched_pattern = pattern
            break

    if not is_core_found:
        return False, "Missing characteristic iridoid/secoiridoid skeleton"

    # Characteristic substituent patterns with more specific connectivity
    key_substituents = [
        ("[C;R2][C;R2]([O;H1])[C;R2]", "hydroxyl on cyclopentane"),
        ("[C;R0][C;R0](=O)[O;H1]", "carboxyl"),
        ("[C;R2][C;H2][O;H1]", "hydroxymethyl"),
        ("[O;R0][C;R0](=O)[C;R0]", "ester"),
        ("[C;R2][C;H3]", "methyl")
    ]
    
    found_groups = []
    for pattern, name in key_substituents:
        pat = Chem.MolFromSmarts(pattern)
        if pat and mol.HasSubstructMatch(pat):
            found_groups.append(name)

    # Ring count check
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 2 or ring_count > 6:
        return False, "Ring count outside typical range for iridoids"

    # Additional checks for secoiridoids
    is_seco = "secoiridoid" if "C;R0][C;R0](=O)" in matched_pattern else "iridoid"
    
    # Check for characteristic oxygen-containing groups
    if len(found_groups) < 2:
        return False, "Insufficient characteristic substituents"

    # Success case
    substituents = ", ".join(found_groups)
    return True, f"Identified as {is_seco} with {substituents}"