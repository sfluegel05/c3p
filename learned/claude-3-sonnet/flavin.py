"""
Classifies: CHEBI:30527 flavin
"""
"""
Classifies: flavin compounds
Based on dimethylisoalloxazine skeleton with N10 substitution
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavin(smiles: str):
    """
    Determines if a molecule is a flavin based on its SMILES string.
    Flavins are derivatives of dimethylisoalloxazine with N10 substitution.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core isoalloxazine pattern - more flexible version
    # Matches the tricyclic system with proper connectivity but allows for different oxidation states
    core_pattern = Chem.MolFromSmarts("""
        [#7]1-,:[#6]2:,-[#6]([#6,#1;!H0]):,-[#6]:,-[#6]:,-[#6]2:,-[#7]([!H0,C])-,:[#6]3:,-[#7]:,-[#6]:,-[#7]:,-[#6]3:,-[#7]1
    """)
    
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Missing isoalloxazine core structure"

    # Check for the characteristic uracil-like ring with carbonyls
    # Allow for different tautomeric forms and oxidation states
    uracil_pattern = Chem.MolFromSmarts("""
        [#7]1-[#6](=[O,S])-[#7]-[#6](=[O,S])-[#6,#7]-[#7]1
    """)
    
    if not mol.HasSubstructMatch(uracil_pattern):
        return False, "Missing characteristic uracil-like ring"

    # Check for N10 substitution (more permissive)
    n10_pattern = Chem.MolFromSmarts("""
        [#6]1:,-[#6]2:,-[#6]:,-[#6]:,-[#6]:,-[#6]2:,-[#7]([#6;!H3])-,:[#6]:,-[#7]:,-[#6]:,-[#7]:,-[#6]:,-[#7]1
    """)
    
    if not mol.HasSubstructMatch(n10_pattern):
        return False, "Missing required N10 substitution"

    # Check for methyl groups on the benzene ring (more flexible)
    # Allow for different positions and modifications
    methyl_pattern = Chem.MolFromSmarts("""
        [#6]1:,-[#6]2:,-[#6]([#6]):,-[#6]([#6]):,-[#6]:,-[#6]2:,-[#7]:,-[#6]:,-[#7]:,-[#6]:,-[#7]:,-[#6]:1
    """)
    
    if not mol.HasSubstructMatch(methyl_pattern):
        return False, "Missing characteristic substitution pattern on benzene ring"

    # Additional validation checks
    
    # Count nitrogens (should have at least 4)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 4:
        return False, "Insufficient number of nitrogen atoms"

    # Check for proper ring system
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 3:
        return False, "Insufficient ring system"

    return True, "Contains isoalloxazine core with proper substitution pattern"

def test_flavins():
    """Test function with known flavins"""
    test_cases = [
        # Riboflavin
        "CC1=C(C)C=C2N(C[C@H](O)[C@H](O)[C@H](O)CO)C3=NC(=O)NC(=O)C3=NC2=C1",
        # FMN
        "C12=NC(NC(C1=NC=3C(N2C[C@@H]([C@@H]([C@@H](COP(=O)(O)O)O)O)O)=CC(=C(C3)C)C)=O)=O",
        # Lumiflavin
        "Cc1cc2nc3c(nc(=O)[nH]c3=O)n(C)c2cc1C",
        # Not a flavin (benzene)
        "c1ccccc1",
        # Not a flavin (purine)
        "n1c2[nH]cnc2nc1"
    ]
    
    for smiles in test_cases:
        result, reason = is_flavin(smiles)
        print(f"SMILES: {smiles}")
        print(f"Is flavin: {result}")
        print(f"Reason: {reason}\n")

if __name__ == "__main__":
    test_flavins()