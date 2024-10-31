from rdkit import Chem
from rdkit.Chem import AllChem

def is_tetracenomycin(smiles: str):
    """
    Determines if a molecule is a tetracenomycin (polyketide based on tetracene ring structure).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a tetracenomycin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for basic ring structure - need 4 fused rings
    rings = mol.GetRingInfo()
    if len(rings.AtomRings()) < 4:
        return False, "Does not contain required minimum of 4 rings"

    # SMARTS pattern for tetracenomycin core structure
    # Four fused 6-membered rings with required oxygenation pattern
    tetracenomycin_pattern = """
        [#6]1[#6][#6]2[#6][#6]3[#6](=[O,#6])[#6]4[#6][#6][#6][#6]4[#6](=[O,#6])[#6]3[#6][#6]2[#6]1
    """
    tetracenomycin_pattern = tetracenomycin_pattern.strip()
    
    pattern_mol = Chem.MolFromSmarts(tetracenomycin_pattern)
    if pattern_mol is None:
        return None, "Invalid SMARTS pattern"
        
    if not mol.HasSubstructMatch(pattern_mol):
        return False, "Does not match tetracenomycin core structure"

    # Check for required oxygenation pattern
    required_patterns = [
        "[OX2H1]", # At least one hydroxyl group
        "[CX3](=O)", # At least one carbonyl
        "[#6]-[OX2H0]-[#6]" # At least one ether
    ]
    
    for pattern in required_patterns:
        pat = Chem.MolFromSmarts(pattern)
        if not mol.HasSubstructMatch(pat):
            return False, f"Missing required functional group: {pattern}"

    # Count oxygens (should have multiple)
    oxygen_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[O]")))
    if oxygen_count < 4:
        return False, "Insufficient oxygen-containing groups"

    return True, "Matches tetracenomycin structure with required functional groups"
# Pr=None
# Recall=0.0