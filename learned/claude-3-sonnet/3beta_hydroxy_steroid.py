"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
"""
Classifies: 3beta-hydroxy steroid
A 3-hydroxy steroid in which the 3-hydroxy substituent is in the beta-position.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 3beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Multiple SMARTS patterns for steroid core to catch different variations
    steroid_patterns = [
        # Basic steroid core (more flexible version)
        "[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~[#6]~4~[#6]~3~[#6]~2~[#6]~1",
        # Alternative pattern with more flexible ring fusion
        "[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~[#6]~4~[#6]~3~[#6]~2~[#6]~1",
        # Pattern allowing for double bonds
        "[#6]1[#6][#6][#6]2[#6][#6][#6]3[#6][#6][#6]4[#6][#6][#6][#6][#6]4[#6]3[#6]2[#6]1"
    ]
    
    has_steroid_core = False
    for pattern in steroid_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_steroid_core = True
            break
            
    if not has_steroid_core:
        return False, "No steroid core structure found"

    # Multiple patterns for 3beta-hydroxy group
    beta_hydroxy_patterns = [
        # Standard 3beta-OH pattern
        '[H][C@@]1[C@@H](O)CC[C@]2',
        # Alternative pattern with different representation
        '[C@@H](O)CC[C@@]1',
        # More general pattern for 3beta-OH
        '[C@@H]1(O)[CH2][CH2]C',
        # Pattern for cyclic systems with 3beta-OH
        '[C@@H](O)[CH2][CH2][C@@]'
    ]
    
    has_beta_hydroxy = False
    for pattern in beta_hydroxy_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_beta_hydroxy = True
            break
            
    if not has_beta_hydroxy:
        return False, "No 3beta-hydroxy group found"

    # Basic structural checks
    # Count rings
    ri = mol.GetRingInfo()
    if ri.NumRings() < 4:
        return False, "Insufficient number of rings for steroid structure"

    # Count carbons (steroids typically have 17+ carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 17:
        return False, "Too few carbons for steroid structure"

    # Count oxygens (should have at least one for the hydroxy group)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 1:
        return False, "No oxygen atoms found"

    # Additional check for sp3 carbons (steroids should have many)
    sp3_carbons = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[C^3]')))
    if sp3_carbons < 10:
        return False, "Too few sp3 carbons for steroid structure"

    return True, "Contains steroid core with 3beta-hydroxy group"