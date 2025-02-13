"""
Classifies: CHEBI:16219 cucurbitacin
"""
"""
Classifies: cucurbitacin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cucurbitacin(smiles: str):
    """
    Determines if a molecule is a cucurbitacin based on its SMILES string.
    Cucurbitacins are tetracyclic triterpenoids with specific structural features.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cucurbitacin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic molecular properties
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 25:
        return False, "Too few carbons for cucurbitacin"
    
    if o_count < 3:
        return False, "Too few oxygens for cucurbitacin"

    # Define core structure SMARTS patterns for cucurbitane skeleton
    core_patterns = [
        # Basic cucurbitane skeleton with various possible arrangements
        "[C]1[C][C]2[C]([C]3[C]([C]4[C]([C][C]3)[C]([C][C]4)(C)C)[C][C]2[C]1",
        # Alternative representation with more specific bond types
        "[CH2]1[CH2][C]2([C]3[C]([C]4[C]([CH2][C]3)[C]([CH2][C]4)(C)C)[C][C]2[CH2]1",
        # Core with common oxygen substitution positions
        "[C]1[C]([OH])[C]2[C]([C]3[C]([C]4[C]([C][C]3)[C]([C][C]4)(C)C)[C][C]2[C]1"
    ]

    # Convert SMARTS to molecules
    core_mols = [Chem.MolFromSmarts(pattern) for pattern in core_patterns]
    core_mols = [mol for mol in core_mols if mol is not None]

    # Check for presence of core structure
    has_core = any(mol.HasSubstructMatch(pattern) for pattern in core_mols)
    if not has_core:
        return False, "Missing cucurbitane core structure"

    # Common functional group patterns in cucurbitacins
    functional_patterns = [
        # C-11 ketone
        "[C]1[C](=O)[C]([C]2[C]1)",
        # C-3 ketone or hydroxyl
        "[C]1[C]([C]2)([OH,=O])",
        # C-16 hydroxyl
        "[C]1[C]([OH])[C]2",
        # C-22 hydroxyl
        "CC(C)(O)C",
        # α,β-unsaturated ketone system
        "C=CC(=O)C",
        # Common side chain pattern
        "CC(C)(O)CC"
    ]

    func_mols = [Chem.MolFromSmarts(pattern) for pattern in functional_patterns]
    func_mols = [mol for mol in func_mols if mol is not None]

    # Count functional group matches
    func_group_count = sum(1 for pattern in func_mols if mol.HasSubstructMatch(pattern))

    # Ring analysis
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient ring count for cucurbitacin"

    # Scoring system
    score = 0
    score += 3 if has_core else 0
    score += func_group_count
    score += 1 if 450 <= rdMolDescriptors.CalcExactMolWt(mol) <= 800 else 0
    score += 1 if o_count >= 5 else 0
    score += 1 if 28 <= c_count <= 35 else 0
    score += 2 if ring_info.NumRings() >= 4 else 0

    # Check for characteristic cucurbitacin features
    if score >= 6 and has_core and func_group_count >= 2:
        return True, "Contains cucurbitane skeleton with characteristic functional groups"
    else:
        return False, f"Insufficient cucurbitacin characteristics (score: {score})"