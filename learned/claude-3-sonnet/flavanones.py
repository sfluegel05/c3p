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

    # Basic flavanone core - chroman-4-one skeleton with aryl at position 2
    # More permissive pattern that allows for various substitutions
    flavanone_core = Chem.MolFromSmarts(
        "[#6]1~2[#6](=[O])[#6][#6]([#6]~[#6]~3~[#6]~[#6]~[#6]~[#6]~[#6]3)O[#6]~2~[#6]~[#6]~[#6]1"
    )
    
    if not mol.HasSubstructMatch(flavanone_core):
        return False, "Missing core flavanone skeleton"

    # Verify saturated C2-C3 bond (distinguishes from flavones)
    # Look specifically for the -O-CH-CH2-C(=O)- pattern
    dihydro_pattern = Chem.MolFromSmarts("O1[CH][CH2]C(=O)")
    if not mol.HasSubstructMatch(dihydro_pattern):
        return False, "Missing required saturated C2-C3 bond"

    # Verify presence of aryl group at C2
    c2_aryl = Chem.MolFromSmarts("O1[CH]([CH2]C(=O))c2ccccc2")
    if not mol.HasSubstructMatch(c2_aryl):
        return False, "Missing required aryl group at C2"

    # Count rings - should have at least 2 (the chroman ring and the aryl ring)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 2:
        return False, "Insufficient ring count"

    # Additional check for ketone at position 4
    ketone_pattern = Chem.MolFromSmarts("[#6]-1-[#6](=O)-[#6]-[#6]-O-[#6]-1")
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "Missing required ketone at position 4"

    return True, "Contains 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton"