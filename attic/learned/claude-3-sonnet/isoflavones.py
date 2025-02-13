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
    
    # 1. Basic isoflavone core - more flexible pattern
    # Matches the basic chromone skeleton with position 3 for aryl attachment
    isoflavone_core = Chem.MolFromSmarts('[#6]1~[#6]~[#6](=O)~c2cccc[c,C]2O1')
    
    # 2. Alternative core pattern to catch variations
    isoflavone_core_alt = Chem.MolFromSmarts('O=C1C=COc2ccccc12')
    
    # 3. Pattern for 3-aryl connection
    aryl_connection = Chem.MolFromSmarts('[#6]1~[#6](~[#6](=O)~c2cccc[c,C]2O1)~[#6]3[#6]~[#6]~[#6]~[#6]~[#6]3')

    if not (mol.HasSubstructMatch(isoflavone_core) or mol.HasSubstructMatch(isoflavone_core_alt)):
        return False, "Missing basic isoflavone core structure"

    # Check ring count (should have at least 3 rings)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 3:
        return False, "Insufficient number of rings"

    # Verify presence of ketone group
    ketone_pattern = Chem.MolFromSmarts('C(=O)')
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "Missing ketone group"

    # Check for allowed atoms
    # Isoflavones typically contain C, H, O and some common substituents
    allowed_atoms = {6, 1, 8, 7, 9, 17, 35, 53}  # C,H,O,N,F,Cl,Br,I
    atom_nums = {atom.GetAtomicNum() for atom in mol.GetAtoms()}
    if not atom_nums.issubset(allowed_atoms):
        return False, "Contains disallowed atoms for isoflavone"

    # Additional checks for common features
    
    # Count oxygen atoms (should have at least 2 for the basic structure)
    o_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8])
    if o_count < 2:
        return False, "Insufficient oxygen atoms for isoflavone structure"

    # Check for aromatic rings (should have at least 2)
    aromatic_rings = 0
    for atom in mol.GetAtoms():
        if atom.GetIsAromatic():
            aromatic_rings += 1
    if aromatic_rings < 8:  # Minimum number of aromatic atoms in basic structure
        return False, "Insufficient aromatic character"

    # If we've passed all checks, this is likely an isoflavone
    return True, "Contains isoflavone core structure with appropriate substitution pattern"