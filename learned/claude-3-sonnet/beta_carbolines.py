"""
Classifies: CHEBI:60834 beta-carbolines
"""
"""
Classifies compounds as beta-carbolines based on their SMILES structure
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_carbolines(smiles: str):
    """
    Determines if a molecule is a beta-carboline based on its SMILES string.
    Beta-carbolines contain a pyrido[3,4-b]indole core structure and their hydrogenated derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-carboline, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core SMARTS patterns for beta-carboline/pyrido[3,4-b]indole structure
    # Pattern 1: Basic pyrido[3,4-b]indole core with flexible saturation
    core_pattern1 = Chem.MolFromSmarts("[#6]1[#6][#6]2[#6]([#6]1)[#7X3]([#6,#1])[#6]3=[#6]2[#7][#6][#6]=[#6]3")
    
    # Pattern 2: Alternative core pattern allowing for spiro fusion
    core_pattern2 = Chem.MolFromSmarts("[#6]1[#6][#6]2[#6]([#6]1)[#7X3][#6]3=[#6]2[#7][#6][#6]=[#6]3")
    
    # Pattern 3: More general pattern for hydrogenated variants
    core_pattern3 = Chem.MolFromSmarts("[#6]1[#6][#6]2[#6]([#6]1)[#7X3][#6]3[#6]2[#7][#6][#6][#6]3")

    # Pattern 4: Specific pattern for common substituted variants seen in examples
    core_pattern4 = Chem.MolFromSmarts("[#6]1[#6][#6]2[#6]([#6]1)[#7X3]([CH3,CH2])[#6]3=[#6]2[#7][#6][#6]=[#6]3")

    # Check for required ring systems
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 3:
        return False, "Insufficient ring systems for beta-carboline structure"

    # Check for any of the core patterns
    if mol.HasSubstructMatch(core_pattern1):
        return True, "Contains pyrido[3,4-b]indole core structure"
    elif mol.HasSubstructMatch(core_pattern2):
        return True, "Contains spiro-fused pyrido[3,4-b]indole core"
    elif mol.HasSubstructMatch(core_pattern3):
        return True, "Contains hydrogenated pyrido[3,4-b]indole core"
    elif mol.HasSubstructMatch(core_pattern4):
        return True, "Contains N-substituted pyrido[3,4-b]indole core"

    # Additional check for tricyclic system with correct atom types
    tricyclic_pattern = Chem.MolFromSmarts("[#6]1[#6][#6]2[#6]([#6]1)[#7][#6]3[#6]2[#7][#6][#6][#6]3")
    if mol.HasSubstructMatch(tricyclic_pattern):
        return True, "Contains beta-carboline-like tricyclic system"

    return False, "Does not contain beta-carboline core structure"