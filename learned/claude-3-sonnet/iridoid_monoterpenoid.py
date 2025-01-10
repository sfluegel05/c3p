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

    # Core iridoid skeleton patterns
    core_patterns = [
        # Classic cyclopentane-pyran fusion
        "[C]1[C][C]2[C]1[O][C][C][C][C]2",
        # Specific iridoid core with oxygen bridge
        "[C]1[C][C]2[O][C][C][C]([C]2[C]1)",
        # Common secoiridoid pattern
        "[C]1[O][C][C]=[C][C]1[C][C][C]",
        # Variation with lactone
        "[C]1[C][C]2[C]1[O][C](=[O])[C][C]2",
        # Specific cyclopentapyran pattern
        "[C]1[C]2[C][C][C]1[O][C][C][C]2",
        # Hemiacetal pattern common in iridoids
        "[C]1[C][C]2[C]1[O][C]([OH])[C][C]2"
    ]

    # Check for core skeleton
    is_core_found = False
    for pattern in core_patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat and mol.HasSubstructMatch(pat):
            is_core_found = True
            break

    # Additional secoiridoid patterns if core not found
    if not is_core_found:
        seco_patterns = [
            # Opened cyclopentane ring with characteristic groups
            "[C][C](=[O])[C][C][O][C]=[C][C]",
            # Common secoiridoid with lactone
            "[C][C](=[O])[O][C][C]=[C][C][C]",
            # Specific pattern with hydroxyl
            "[C][C]([OH])[C][C][O][C]=[C][C]"
        ]
        for pattern in seco_patterns:
            pat = Chem.MolFromSmarts(pattern)
            if pat and mol.HasSubstructMatch(pat):
                is_core_found = True
                break

    if not is_core_found:
        return False, "Missing characteristic iridoid/secoiridoid skeleton"

    # Check for characteristic substituents in specific positions
    key_substituents = [
        ("[C][C]([OH])[C]", "hydroxyl on cyclopentane"),
        ("[C][C](=[O])[OH]", "carboxyl"),
        ("[C][CH2][OH]", "hydroxymethyl"),
        ("[O][C](=[O])[C]", "ester"),
        ("[C][CH3]", "methyl")
    ]
    
    found_groups = []
    for pattern, name in key_substituents:
        pat = Chem.MolFromSmarts(pattern)
        if pat and mol.HasSubstructMatch(pat):
            found_groups.append(name)

    if len(found_groups) < 2:
        return False, "Insufficient characteristic substituents"

    # Exclude compounds that are primarily sugars
    sugar_patterns = [
        # Pyranose sugar pattern
        "[OH]C1[CH]([OH])[CH]([OH])[CH]([OH])[CH]([OH])O1",
        # Furanose sugar pattern
        "[OH]C1[CH]([OH])[CH]([OH])[CH]([OH])O1",
        # Multiple connected sugars
        "O1[CH]([CH2][OH])[CH]([OH])[CH]([OH])[CH]([OH])[CH]1O[CH]2O[CH]([CH2][OH])[CH]([OH])[CH]([OH])[CH]2[OH]"
    ]
    
    for pattern in sugar_patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat and mol.HasSubstructMatch(pat):
            # Check if the sugar part is dominant
            if o_count > 6 and rdMolDescriptors.CalcNumRings(mol) < 3:
                return False, "Structure appears to be primarily a glycoside"

    # Ring count check
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 2 or ring_count > 6:
        return False, "Ring count outside typical range for iridoids"

    # Success case
    structure_type = "secoiridoid" if any("seco" in p for p in found_groups) else "iridoid"
    substituents = ", ".join(found_groups)
    return True, f"Identified as {structure_type} with {substituents}"