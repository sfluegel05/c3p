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

    # Basic checks
    if mol.GetNumAtoms() < 15:  # Minimum atoms for isoflavone core
        return False, "Too few atoms for isoflavone"

    # Core structure patterns - break it down into parts
    
    # 1. Chromone core (benzopyran-4-one)
    # Note: Using simpler pattern that matches both aromatic and non-aromatic forms
    chromone_pattern = Chem.MolFromSmarts('O=C1c2ccccc2OC=C1')
    if chromone_pattern is None:
        return None, "Invalid SMARTS pattern for chromone"
    if not mol.HasSubstructMatch(chromone_pattern):
        return False, "Missing chromone core structure"

    # 2. Check for aryl group at position 3
    # This pattern looks for the connection between chromone and phenyl ring
    isoflavone_pattern = Chem.MolFromSmarts('O=C1c2ccccc2OC(c3ccccc3)=C1')
    if isoflavone_pattern is None:
        return None, "Invalid SMARTS pattern for isoflavone"
    if not mol.HasSubstructMatch(isoflavone_pattern):
        return False, "Missing aryl group at position 3 (not an isoflavone)"

    # 3. Ring count check (should have at least 3 rings)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 3:
        return False, "Insufficient number of rings"

    # Additional validation - check atom types
    # Isoflavones should only contain C, H, O (and possibly other common substituents)
    allowed_atoms = {6, 1, 8, 7, 9, 17, 35, 53}  # C,H,O,N,F,Cl,Br,I
    atom_nums = {atom.GetAtomicNum() for atom in mol.GetAtoms()}
    if not atom_nums.issubset(allowed_atoms):
        return False, "Contains disallowed atoms for isoflavone"

    # Success - all patterns matched
    return True, "Contains 3-aryl-1-benzopyran-4-one skeleton with appropriate substitution pattern"