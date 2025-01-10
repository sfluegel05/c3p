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
    Isoflavones have a 3-aryl-1-benzopyran-4-one (3-phenylchromen-4-one) skeleton.

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

    # Basic checks
    if mol.GetNumAtoms() < 15:  # Minimum atoms for isoflavone core
        return False, "Too few atoms for isoflavone"

    # Core structure patterns
    
    # 1. Basic isoflavone core with aromatic or non-aromatic forms
    # Matches the core structure including:
    # - benzopyran system (fused rings)
    # - ketone at position 4
    # - connection point for aryl group at position 3
    isoflavone_core = Chem.MolFromSmarts('[#6]1=[#6]-[#6](=O)-c2[c,C]([c,C][c,C][c,C][c,C]2)O1')
    if not mol.HasSubstructMatch(isoflavone_core):
        return False, "Missing basic isoflavone core structure"

    # 2. Check for aryl group at position 3
    # More flexible pattern that allows for substituted phenyl rings
    aryl_pattern = Chem.MolFromSmarts('[#6]1=[#6]-[#6](=O)-c2[c,C]([c,C][c,C][c,C][c,C]2)O1-[#6]3=[#6][#6]=[#6][#6]=[#6]3')
    if not mol.HasSubstructMatch(aryl_pattern):
        return False, "Missing or incorrect aryl group at position 3"

    # 3. Ring count check (should have at least 3 rings - two from chromone core and one from aryl group)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 3:
        return False, "Insufficient number of rings"

    # 4. Check for proper connectivity
    # Ensure the aryl group is connected at position 3 of the chromone
    connectivity_pattern = Chem.MolFromSmarts('[#6]1=COc2[c,C][c,C][c,C][c,C]c2[#6](=O)1')
    if not mol.HasSubstructMatch(connectivity_pattern):
        return False, "Incorrect connectivity pattern"

    # 5. Check for allowed atoms
    # Isoflavones typically contain C, H, O and some common substituents
    allowed_atoms = {6, 1, 8, 7, 9, 17, 35, 53}  # C,H,O,N,F,Cl,Br,I
    atom_nums = {atom.GetAtomicNum() for atom in mol.GetAtoms()}
    if not atom_nums.issubset(allowed_atoms):
        return False, "Contains disallowed atoms for isoflavone"

    # Additional structural validation
    # Check for presence of ketone group
    ketone_pattern = Chem.MolFromSmarts('C(=O)')
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "Missing ketone group"

    return True, "Contains 3-aryl-1-benzopyran-4-one skeleton with appropriate substitution pattern"