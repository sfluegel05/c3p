"""
Classifies: CHEBI:38757 isoflavones
"""
"""
Classifies: CHEBI:47919 isoflavones
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_isoflavones(smiles: str):
    """
    Determines if a molecule is an isoflavone based on its SMILES string.
    Isoflavones have a 3-aryl-1-benzopyran-4-one skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core structure of isoflavone (3-aryl-1-benzopyran-4-one)
    # [#6]:benzene ring connected to oxygen in chromone
    # [#6]=O: carbonyl group
    # [#6]1: carbon forming the pyrone ring
    # [#6]: carbon at position 3 where aryl group attaches
    isoflavone_core = Chem.MolFromSmarts('[#6]2~[#6]:,[#6]~[#6]:,[#6]~[#6]:,[#6]~[#6]2~[#8]~[#6]1~[#6](=[#8])~[#6](~[#6]3:[#6]:[#6]:[#6]:[#6]:[#6]:3)~[#6]~[#8]~1')
    
    if not mol.HasSubstructMatch(isoflavone_core):
        return False, "Missing isoflavone core structure (3-aryl-1-benzopyran-4-one)"

    # Verify presence of ketone group in correct position
    ketone_pattern = Chem.MolFromSmarts('O=C1C=COC2=CC=CC=C12')
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "Missing ketone group in correct position"

    # Check for aryl group at position 3
    aryl_at_3_pattern = Chem.MolFromSmarts('O=C1C=C(c2ccccc2)Oc2ccccc12')
    if not mol.HasSubstructMatch(aryl_at_3_pattern):
        return False, "Missing aryl group at position 3"

    # Count rings to ensure proper structure
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 3:  # Must have at least 3 rings (2 from chromone + 1 from aryl)
        return False, "Insufficient number of rings"

    # Additional check to ensure proper fusion of rings
    fused_rings_pattern = Chem.MolFromSmarts('c2cc1occ(c3ccccc3)c(=O)c1cc2')
    if not mol.HasSubstructMatch(fused_rings_pattern):
        return False, "Incorrect ring fusion pattern"

    return True, "Contains 3-aryl-1-benzopyran-4-one skeleton with proper substitution pattern"