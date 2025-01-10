"""
Classifies: CHEBI:28863 flavanones
"""
"""
Classifies: CHEBI:24064 flavanone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavanones(smiles: str):
    """
    Determines if a molecule is a flavanone based on its SMILES string.
    Flavanones have a 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavanone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core pattern 1: Basic chromanone skeleton with flexible aryl attachment
    # Allows for substitutions and variations in aromatic systems
    core_pattern = Chem.MolFromSmarts(
        "[#6]~1~2~[#6](=[O])~[#6]~[#6]~[#6](~[#6]~3~[#6,#7]~[#6,#7]~[#6,#7]~[#6,#7]~[#6,#7]3)~O~[#6]~2~[#6]~[#6]~[#6]~1"
    )
    
    # Core pattern 2: Alternative representation focusing on the key bonds
    alt_core = Chem.MolFromSmarts(
        "O1[CH][CH2]C(=O)c2c1cccc2"
    )
    
    if not (mol.HasSubstructMatch(core_pattern) or mol.HasSubstructMatch(alt_core)):
        return False, "Missing core flavanone skeleton"

    # Check for ketone at position 4 (more permissive pattern)
    ketone_pattern = Chem.MolFromSmarts("[#6]-1-[#6](=O)-[#6]-[#6]-O-[#6]-1")
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "Missing required ketone at position 4"

    # Check for saturated C2-C3 bond (more permissive pattern)
    # Allows for substitutions at these positions
    dihydro_pattern = Chem.MolFromSmarts("O1[#6][#6]C(=O)")
    if not mol.HasSubstructMatch(dihydro_pattern):
        return False, "Missing required saturated C2-C3 bond"

    # Verify presence of any aromatic ring system at C2 position
    # More permissive pattern that allows for substituted aryl groups
    c2_aryl = Chem.MolFromSmarts("O1[#6]([#6]C(=O))~[#6]~2~[#6,#7]~[#6,#7]~[#6,#7]~[#6,#7]~[#6,#7]2")
    if not mol.HasSubstructMatch(c2_aryl):
        return False, "Missing required aryl group at C2"

    # Count rings - should have at least 2 (the chroman ring and the aryl ring)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 2:
        return False, "Insufficient ring count"

    # Additional structural check for the basic flavanone framework
    framework = Chem.MolFromSmarts("O1CC(=O)c2ccccc2[#6]1")
    if not mol.HasSubstructMatch(framework):
        return False, "Invalid flavanone framework"

    return True, "Contains 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton"