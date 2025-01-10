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

    # Basic beta-carboline core pattern - more flexible version
    # Matches the essential tricyclic system with correct connectivity
    core_pattern = Chem.MolFromSmarts("[#6,#7]1[#6,#7]2[#6]3[#6][#6][#6][#6]3[#7][#6]2=[#6][#6][#6]1")
    
    # Alternative pattern for hydrogenated variants
    hydro_pattern = Chem.MolFromSmarts("[#6,#7]1[#6,#7]2[#6]3[#6][#6][#6][#6]3[#7][#6]2[#6][#6][#6]1")
    
    # Pattern for spiro-fused variants commonly seen in examples
    spiro_pattern = Chem.MolFromSmarts("[#6,#7]1[#6,#7]2[#6]3[#6][#6][#6][#6]3[#7][#6]2[#6]4([#6][#6]1)[#6,#7][#6,#7][#6,#7]4")

    # Check for required ring systems
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 3:
        return False, "Insufficient ring systems for beta-carboline structure"

    # Count nitrogen atoms in molecule
    n_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7])
    if n_count < 2:
        return False, "Insufficient nitrogen atoms for beta-carboline structure"

    # Check for any of the core patterns
    if mol.HasSubstructMatch(core_pattern):
        return True, "Contains beta-carboline core structure"
    elif mol.HasSubstructMatch(hydro_pattern):
        return True, "Contains hydrogenated beta-carboline core"
    elif mol.HasSubstructMatch(spiro_pattern):
        return True, "Contains spiro-fused beta-carboline core"

    # Additional check for common substitution pattern
    substituted_pattern = Chem.MolFromSmarts("[#6,#7]1[#6,#7]2[#6]3[#6][#6]([O,N,C])[#6][#6]3[#7]([C,H])[#6]2[#6][#6][#6]1")
    if mol.HasSubstructMatch(substituted_pattern):
        return True, "Contains substituted beta-carboline core"

    return False, "Does not contain beta-carboline core structure"