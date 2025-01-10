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

    # Define steroid core patterns - try multiple patterns to catch different representations
    steroid_patterns = [
        # Basic steroid core with flexible bonds
        "[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~[#6]~[#6]~3~[#6]~2~[#6]~[#6]~1",
        # Alternative pattern focusing on the ABCD ring connectivity
        "[#6]1~[#6]~[#6]~[#6]2~[#6]~[#6]~[#6]3~[#6]~[#6]~[#6]4~[#6]~[#6]~[#6](~[#6]~4)~[#6]~3~[#6]~2~1"
    ]
    
    has_steroid_core = False
    for pattern in steroid_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_steroid_core = True
            break
            
    if not has_steroid_core:
        return False, "No steroid core structure found"

    # Look for 11β-hydroxy group with correct stereochemistry
    # This pattern specifically looks for the beta-oriented OH at position 11
    # The pattern includes the key carbon atoms around position 11 to ensure correct connectivity
    beta_hydroxy_patterns = [
        # Pattern for 11β-OH with explicit stereochemistry
        "[#6]~1~[#6]~[#6]~[#6]~2~[#6](@[H])~[#6]~[C@@]([H])(O)~[#6]~[#6]~2~[#6]~1",
        # Alternative pattern for cases where stereochemistry is implicit
        "[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[C@H](O)~[#6]~[#6]~2~[#6]~1"
    ]
    
    has_beta_hydroxy = False
    for pattern in beta_hydroxy_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_beta_hydroxy = True
            break
            
    if not has_beta_hydroxy:
        return False, "No 11-beta hydroxy group found with correct stereochemistry"

    # Basic validation of structure
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings for steroid structure"

    # Count carbons and oxygens
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    num_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if num_carbons < 19:
        return False, f"Too few carbons ({num_carbons}) for steroid structure"
    
    if num_oxygens < 1:
        return False, "Must have at least one oxygen atom"

    # Look for hydroxy groups
    hydroxy_pattern = Chem.MolFromSmarts("[OH]")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No hydroxy groups found"

    return True, "Contains steroid core with 11-beta hydroxy group"