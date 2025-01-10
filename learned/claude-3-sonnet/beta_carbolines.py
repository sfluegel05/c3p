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

    # Basic beta-carboline core pattern (aromatic)
    # Matches the tricyclic system with correct nitrogen positions
    core_pattern = Chem.MolFromSmarts("c1ccc2c(c1)[nH]c3c2c(=*)nc[c,n]3")
    
    # Hydrogenated beta-carboline core pattern
    # Allows for reduced forms of the core
    hydro_pattern = Chem.MolFromSmarts("C1CC=C2C(=C1)NC3=C2C(CNC3)=*")
    
    # Pattern for N-substituted variants
    n_sub_pattern = Chem.MolFromSmarts("c1ccc2c(c1)n(C)c3c2c(C*)ncc3")
    
    # Pattern for spiro variants at position 4
    spiro_pattern = Chem.MolFromSmarts("c1ccc2c(c1)[nH]c3c2C4(CNC3)CC4")
    
    # Check ring count
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 3:
        return False, "Insufficient ring systems for beta-carboline structure"

    # Count nitrogen atoms
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 2:
        return False, "Insufficient nitrogen atoms for beta-carboline structure"

    # Check for any of the core patterns
    if mol.HasSubstructMatch(core_pattern):
        return True, "Contains beta-carboline core structure"
    
    if mol.HasSubstructMatch(hydro_pattern):
        return True, "Contains hydrogenated beta-carboline core"
    
    if mol.HasSubstructMatch(n_sub_pattern):
        return True, "Contains N-substituted beta-carboline core"
    
    if mol.HasSubstructMatch(spiro_pattern):
        return True, "Contains spiro-fused beta-carboline core"

    # Additional check for 7-substituted variants (common in examples)
    sub_pattern = Chem.MolFromSmarts("c1cc(*)c2c(c1)nc3c2c(C*)ncc3")
    if mol.HasSubstructMatch(sub_pattern):
        return True, "Contains substituted beta-carboline core"

    return False, "Does not contain beta-carboline core structure"