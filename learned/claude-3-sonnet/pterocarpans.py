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

    # Core pterocarpan skeleton SMARTS pattern
    # This pattern specifically matches the 6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene skeleton
    # [#6@H] represents the sp3 carbons with specific stereochemistry
    pterocarpan_core = Chem.MolFromSmarts(
        '[#6]1-2=[#6]-[#6](=[#6]-[#6]=[#6]-1)-[#8]-[#6@]3-[#6]-[#8]-[#6]4=[#6]-[#6]=[#6]-[#6]=[#6]-4-[#6]-2-3'
    )
    
    # Alternative pattern that's more flexible for substituted derivatives
    pterocarpan_core_alt = Chem.MolFromSmarts(
        '[#6]1-2~[#6]~[#6](~[#6]~[#6]~[#6]~1)~[#8]~[#6]3~[#6]~[#8]~[#6]4~[#6]~[#6]~[#6]~[#6]~[#6]~4~[#6]~2~3'
    )

    if not (mol.HasSubstructMatch(pterocarpan_core) or mol.HasSubstructMatch(pterocarpan_core_alt)):
        return False, "Missing characteristic pterocarpan core structure"

    # Check for the specific fusion pattern with two sp3 carbons
    fusion_pattern = Chem.MolFromSmarts('[#6X4]1-[#8]-[#6]-c2ccccc2-[#6X4]1-[#8]')
    if not mol.HasSubstructMatch(fusion_pattern):
        return False, "Missing required sp3 carbons at fusion points"

    # Verify the presence of both oxygen bridges in correct positions
    oxygen_bridge_pattern = Chem.MolFromSmarts('[#6]1-[#8]-[#6]-[#6]-[#8]-[#6]-1')
    if not mol.HasSubstructMatch(oxygen_bridge_pattern):
        return False, "Missing characteristic oxygen bridge pattern"

    # Check for the benzofuran system
    benzofuran_pattern = Chem.MolFromSmarts('c1cccc2c1OC-[#6]2')
    if not mol.HasSubstructMatch(benzofuran_pattern):
        return False, "Missing benzofuran system"

    # Check for the chromane system
    chromane_pattern = Chem.MolFromSmarts('c1cccc2c1OC[#6]2')
    if not mol.HasSubstructMatch(chromane_pattern):
        return False, "Missing chromane system"

    # Additional check to exclude coumestans (which would have a C=O group)
    carbonyl_pattern = Chem.MolFromSmarts('[#6]1=[#8]-[#6]2=,:[#6]-[#6]=,:[#6]-[#6]=,:[#6]-2-[#8]-1')
    if mol.HasSubstructMatch(carbonyl_pattern):
        return False, "Structure appears to be a coumestan rather than pterocarpan"

    return True, "Contains characteristic 6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene skeleton"