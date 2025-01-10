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

    # Molecular weight check (typical range for iridoids)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150 or mol_wt > 800:
        return False, "Molecular weight outside typical range for iridoids"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 8 or c_count > 30:
        return False, "Carbon count outside typical range for iridoids"
    
    if o_count < 2:
        return False, "Insufficient oxygen atoms for iridoid structure"

    # Define various iridoid skeleton patterns
    patterns = [
        # Classic iridoid skeleton
        "[C]1[C][C][C]2[C]1[O][C][C][C][C]2",
        # Seco-iridoid skeleton
        "[C]-[C]-[C]-[C]1[O][C][C][C][C]1",
        # Modified iridoid skeleton with oxygen bridge
        "[C]1[C][C]2[O][C][C][C][C]2[O]1",
        # Rearranged iridoid skeleton
        "[C]1[C]2[C][C][C]1[O][C][C][C]2",
        # Cyclopentane fused to pyran (more general)
        "[C]1[C][C]2[C][C]1[O][C][C][C]2"
    ]

    is_iridoid_skeleton = False
    for pattern in patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat and mol.HasSubstructMatch(pat):
            is_iridoid_skeleton = True
            break

    if not is_iridoid_skeleton:
        # Check for secoiridoid patterns
        seco_patterns = [
            # Common secoiridoid patterns
            "[C]-[C]-[C]-[C]-[C]-[O]-[C]",
            "[C]-[C](=O)-[C]-[C]-[O]-[C]",
            "[C]-[C](=O)-[C]-[C]-[C]-[O]"
        ]
        
        for pattern in seco_patterns:
            pat = Chem.MolFromSmarts(pattern)
            if pat and mol.HasSubstructMatch(pat):
                is_iridoid_skeleton = True
                break

    if not is_iridoid_skeleton:
        return False, "Missing characteristic iridoid/secoiridoid skeleton"

    # Check for characteristic substituent patterns
    substituent_patterns = [
        (Chem.MolFromSmarts("[CH2][OH]"), "hydroxymethyl"),
        (Chem.MolFromSmarts("[C](=O)[OH]"), "carboxyl"),
        (Chem.MolFromSmarts("[OH]"), "hydroxyl"),
        (Chem.MolFromSmarts("[C](=O)[O][C]"), "ester"),
        (Chem.MolFromSmarts("[CH3]"), "methyl")
    ]
    
    found_groups = []
    for pattern, name in substituent_patterns:
        if pattern and mol.HasSubstructMatch(pattern):
            found_groups.append(name)
    
    if len(found_groups) < 2:
        return False, "Insufficient characteristic substituents"

    # Ring count check
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 2 or ring_count > 6:
        return False, "Ring count outside typical range for iridoids"

    # Exclude compounds that are likely to be simple glycosides
    sugar_pattern = Chem.MolFromSmarts("[OH]C1[CH]([OH])[CH]([OH])[CH]([OH])[CH]([OH])[CH]1[OH]")
    if sugar_pattern and mol.HasSubstructMatch(sugar_pattern):
        # Additional check to ensure it's not just a simple glycoside
        if ring_count == 2 and o_count > 5:
            return False, "Likely a simple glycoside"

    # Success case
    structure_type = "secoiridoid" if not any("cyclopentane" in p for p in patterns) else "iridoid"
    substituents = ", ".join(found_groups)
    return True, f"Identified as {structure_type} with {substituents} groups"