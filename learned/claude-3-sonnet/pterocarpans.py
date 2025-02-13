"""
Classifies: CHEBI:26377 pterocarpans
"""
"""
Classifies: pterocarpans
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_pterocarpans(smiles: str):
    """
    Determines if a molecule is a pterocarpan based on its SMILES string.
    Pterocarpans have a characteristic 6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pterocarpan, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic pterocarpan core skeleton - more flexible pattern
    # Represents the basic connectivity without being too restrictive
    core_pattern = Chem.MolFromSmarts(
        '[#6]1~2~[#6]~[#6]~[#6]~[#6]~[#6]~1~[#8]~[#6]~3~[#6]~[#8]~[#6]~4~[#6]~[#6]~[#6]~[#6]~[#6]~4~[#6]~2~3'
    )
    
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Missing basic pterocarpan skeleton"

    # Check for the two sp3 carbons at the junction (6a,11a positions)
    # Using a more specific pattern for the fusion points
    junction_pattern = Chem.MolFromSmarts('[#6X4]1-[#8]-[#6]-[#6]-2-[#6]-[#6]-[#6]-[#6]-[#6]-2-[#6X4]-1')
    if not mol.HasSubstructMatch(junction_pattern):
        return False, "Missing required sp3 carbons at 6a,11a positions"

    # Check for the characteristic oxygen bridges
    # More flexible pattern that allows for various substitutions
    oxygen_bridges = Chem.MolFromSmarts('[#6]1-[#8]-[#6]-[#6]2-[#8]-[#6]-[#6]-[#6]-1-2')
    if not mol.HasSubstructMatch(oxygen_bridges):
        return False, "Missing characteristic oxygen bridge pattern"

    # Exclude molecules with carbonyl groups in the core structure (coumestans)
    carbonyl_pattern = Chem.MolFromSmarts('[#6]1=[#8]-[#6]2=,:c-[#6]=,:[#6]-[#6]=,:[#6]-2-[#8]-1')
    if mol.HasSubstructMatch(carbonyl_pattern):
        return False, "Contains carbonyl group characteristic of coumestans"

    # Additional check for proper ring connectivity
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings for pterocarpan structure"

    # Check for proper aromatic systems
    aromatic_pattern = Chem.MolFromSmarts('c1cccc2c1-[#6]-[#8]-[#6]-2')
    if len(mol.GetSubstructMatches(aromatic_pattern)) < 2:
        return False, "Missing required aromatic systems"

    return True, "Contains characteristic pterocarpan skeleton with proper ring fusion and oxygen bridges"