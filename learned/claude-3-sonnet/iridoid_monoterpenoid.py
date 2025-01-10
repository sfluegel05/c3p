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
    Iridoids typically have a cyclopentane ring fused to a six-membered oxygen heterocycle.

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

    # Basic molecular properties check
    num_atoms = mol.GetNumAtoms()
    if num_atoms < 10 or num_atoms > 50:
        return False, "Molecule size outside typical range for iridoid monoterpenoids"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 8:
        return False, "Too few carbons for an iridoid monoterpenoid"
    
    if o_count < 1:
        return False, "No oxygen atoms found"

    # Look for cyclopentane ring (5-membered carbon ring)
    cyclopentane_pattern = Chem.MolFromSmarts("[C]1[C][C][C][C]1")
    
    # Look for pyran ring (6-membered ring with oxygen)
    pyran_pattern = Chem.MolFromSmarts("[O]1[C][C][C][C][C]1")
    
    # Look for broken cyclopentane (for secoiridoids)
    seco_pattern = Chem.MolFromSmarts("[C]-[C]-[C]-[C]-[C]")
    
    # Check for presence of either intact cyclopentane or seco form
    has_cyclopentane = mol.HasSubstructMatch(cyclopentane_pattern)
    has_pyran = mol.HasSubstructMatch(pyran_pattern)
    has_seco = mol.HasSubstructMatch(seco_pattern)
    
    if not (has_cyclopentane or has_seco):
        return False, "Missing required cyclopentane or seco-cyclopentane structure"
    
    if not has_pyran:
        # Look for modified oxygen-containing heterocycle
        modified_pyran = Chem.MolFromSmarts("[O]1[C][C,O][C,O][C][C]1")
        if not mol.HasSubstructMatch(modified_pyran):
            return False, "Missing required pyran or modified oxygen heterocycle"

    # Check for typical substituents often found in iridoids
    common_groups = [
        (Chem.MolFromSmarts("[OH]"), "hydroxyl"),
        (Chem.MolFromSmarts("[C](=O)[OH]"), "carboxyl"),
        (Chem.MolFromSmarts("[CH3]"), "methyl"),
        (Chem.MolFromSmarts("[CH2][OH]"), "hydroxymethyl")
    ]
    
    found_groups = []
    for pattern, name in common_groups:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            found_groups.append(name)
    
    if not found_groups:
        return False, "Missing typical iridoid substituents"

    # Calculate ring count
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 2:
        return False, "Insufficient ring count for iridoid structure"
    
    # Success case
    structure_type = "seco-iridoid" if has_seco else "iridoid"
    substituents = ", ".join(found_groups)
    return True, f"Identified as {structure_type} with {substituents} groups"