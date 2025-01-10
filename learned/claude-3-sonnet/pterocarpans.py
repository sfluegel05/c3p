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

    # Basic pterocarpan core skeleton
    # More permissive pattern focusing on connectivity
    core_pattern = Chem.MolFromSmarts(
        '[#6]1~2~[#6]~[#6]~[#6]~[#6]~[#6]~1~[#8]~[#6]~3~[#6]~[#8]~[#6]([#6]~2~3])'
    )
    
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Missing basic pterocarpan skeleton"

    # Check for the characteristic benzofuran system
    benzofuran = Chem.MolFromSmarts('c1cccc2c1OCC2')
    if not mol.HasSubstructMatch(benzofuran):
        return False, "Missing benzofuran system"

    # Check for the characteristic chromene system
    chromene = Chem.MolFromSmarts('c1cccc2c1OC-C2')
    if not mol.HasSubstructMatch(chromene):
        return False, "Missing chromene system"

    # Verify the presence of two sp3 carbons at the junction (6a,11a positions)
    junction = Chem.MolFromSmarts('[CH1X4]1-[CH2X4]-O-c2ccccc2-[CH1X4]1')
    if not mol.HasSubstructMatch(junction):
        # Try alternative pattern for substituted junctions
        junction_alt = Chem.MolFromSmarts('[CX4]1-[CH2X4]-O-c2ccccc2-[CX4]1')
        if not mol.HasSubstructMatch(junction_alt):
            return False, "Missing required sp3 carbons at 6a,11a positions"

    # Check ring count (should have at least 4 rings)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings"

    # Exclude molecules that are actually coumestans (oxidized pterocarpans)
    coumestan = Chem.MolFromSmarts('O=C1Oc2ccccc2-c2ccc3ccccc3o12')
    if mol.HasSubstructMatch(coumestan):
        return False, "Structure is a coumestan (oxidized pterocarpan)"

    # Additional check for proper fusion pattern
    fusion = Chem.MolFromSmarts('c1ccc2c(c1)OC[C@H]1[C@H]2Oc2ccccc21')
    fusion_alt = Chem.MolFromSmarts('c1ccc2c(c1)OC[C@@H]1[C@@H]2Oc2ccccc21')
    
    if not (mol.HasSubstructMatch(fusion) or mol.HasSubstructMatch(fusion_alt)):
        # Try a more general fusion pattern without stereochemistry
        fusion_general = Chem.MolFromSmarts('c1ccc2c(c1)OCC1C2Oc2ccccc21')
        if not mol.HasSubstructMatch(fusion_general):
            return False, "Incorrect ring fusion pattern"

    return True, "Contains characteristic pterocarpan skeleton with proper ring fusion and oxygen bridges"