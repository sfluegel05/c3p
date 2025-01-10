"""
Classifies: CHEBI:35346 11beta-hydroxy steroid
"""
"""
Classifies: 11beta-hydroxy steroids
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_11beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is an 11beta-hydroxy steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an 11beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic steroid core pattern that's more flexible
    # Allows for both single and double bonds in the skeleton
    steroid_core = Chem.MolFromSmarts("[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~[#6]~[#6]~3~[#6]~2~[#6]~[#6]~1")
    
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # More specific pattern for 11β-hydroxy group
    # This pattern looks for the characteristic 11β-OH configuration in the context of rings C and D
    hydroxy_11beta_pattern = Chem.MolFromSmarts("[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[C@@H](O)~[#6]~[#6]~2~[#6]~1")
    
    if not mol.HasSubstructMatch(hydroxy_11beta_pattern):
        return False, "No 11-beta hydroxy group found"

    # Validate ring structure
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings for steroid structure"

    # Count carbons and oxygens
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    num_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if num_carbons < 19:
        return False, f"Too few carbons ({num_carbons}) for steroid structure"
    
    if num_oxygens < 2:
        return False, "Insufficient oxygen atoms for 11-beta-hydroxy steroid"

    # Check for common steroid features
    ketone_pattern = Chem.MolFromSmarts("C(=O)")
    has_ketone = mol.HasSubstructMatch(ketone_pattern)
    
    hydroxy_pattern = Chem.MolFromSmarts("[#6][OH]")
    num_hydroxy = len(mol.GetSubstructMatches(hydroxy_pattern))
    
    if num_hydroxy < 1:
        return False, "No hydroxy groups found"

    # Additional check for common substituents at C-17 and C-21 positions
    c17_c21_pattern = Chem.MolFromSmarts("[#6]-[#6](=O)-[#6,#8]")
    if not mol.HasSubstructMatch(c17_c21_pattern) and not has_ketone:
        return False, "Missing typical steroid oxygenated substituents"

    # Check molecular weight
    mol_weight = Chem.Descriptors.ExactMolWt(mol)
    if mol_weight < 250 or mol_weight > 800:  # Increased upper limit to account for larger derivatives
        return False, f"Molecular weight {mol_weight:.1f} outside typical range for steroids"

    # If all checks pass, it's likely an 11β-hydroxy steroid
    return True, "Contains steroid core with 11-beta hydroxy group and appropriate substituents"