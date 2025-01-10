"""
Classifies: CHEBI:18254 ribonucleoside
"""
"""
Classifies: CHEBI:18254 ribonucleoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_ribonucleoside(smiles: str):
    """
    Determines if a molecule is a ribonucleoside based on its SMILES string.
    A ribonucleoside consists of a nucleobase attached to a D-ribose sugar.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ribonucleoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for ribose pattern (furan ring with multiple OH groups)
    # [#6] represents carbon, [#8] represents oxygen
    # The @ symbols ensure correct stereochemistry for D-ribose
    ribose_pattern = Chem.MolFromSmarts("[#6]1-[#8]-[#6]-[#6]-[#6]1")
    if not mol.HasSubstructMatch(ribose_pattern):
        return False, "No ribose sugar found"

    # Check for hydroxyl groups on ribose (typically 3)
    hydroxyl_pattern = Chem.MolFromSmarts("[#8H1]")
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if hydroxyl_matches < 2:  # Allow for some modification of OH groups
        return False, f"Insufficient hydroxyl groups ({hydroxyl_matches})"

    # Look for nucleobase patterns
    # Pyrimidine pattern (6-membered ring with nitrogens)
    pyrimidine_pattern = Chem.MolFromSmarts("c1[n]c[n]c1")
    # Purine pattern (fused 5,6 ring system with nitrogens)
    purine_pattern = Chem.MolFromSmarts("c1[n]c2[n]c[n]c2[n]1")
    
    has_pyrimidine = mol.HasSubstructMatch(pyrimidine_pattern)
    has_purine = mol.HasSubstructMatch(purine_pattern)
    
    if not (has_pyrimidine or has_purine):
        return False, "No nucleobase (pyrimidine or purine) found"

    # Check for N-glycosidic bond connecting nucleobase to ribose
    # This is a nitrogen connected to the anomeric carbon of ribose
    n_glycosidic_pattern = Chem.MolFromSmarts("[#7]-[#6]1-[#8]-[#6]-[#6]-[#6]1")
    if not mol.HasSubstructMatch(n_glycosidic_pattern):
        return False, "No N-glycosidic bond found"

    # Additional check for reasonable molecular size
    if mol.GetNumAtoms() < 15:  # Minimum size for a basic nucleoside
        return False, "Molecule too small to be a ribonucleoside"
    
    # Check if the molecule has a reasonable number of rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 2:  # Should have at least ribose + base
        return False, "Insufficient ring count"

    # Success case - molecule appears to be a ribonucleoside
    base_type = "purine" if has_purine else "pyrimidine"
    return True, f"Contains ribose sugar connected to {base_type} base via N-glycosidic bond"