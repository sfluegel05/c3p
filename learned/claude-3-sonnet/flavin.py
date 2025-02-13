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

    # Basic tricyclic system pattern (pteridine fused with benzene)
    # More permissive pattern that captures the essential connectivity
    core_pattern = Chem.MolFromSmarts("[#7]1~[#6]2~[#6]~[#6]~[#6]~[#6]2~[#7]~[#6]3~[#7]~[#6]~[#7]~[#6]3~[#7]1")
    if core_pattern is None:
        return None, "Error in SMARTS pattern"
    
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Missing core tricyclic flavin system"

    # Check for the two carbonyl groups in the uracil-like ring
    carbonyl_pattern = Chem.MolFromSmarts("[#7]1~[#6](=[O])~[#7]~[#6](=[O])~[#6]~[#7]1")
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return False, "Missing characteristic carbonyl groups"

    # Check for N10 substitution pattern
    n10_pattern = Chem.MolFromSmarts("[#6]1~[#6]2~[#6]~[#6]~[#6]~[#6]2~[#7]([#6,#1;!H0])~[#6]3~[#7]~[#6]~[#7]~[#6]3~[#7]1")
    if not mol.HasSubstructMatch(n10_pattern):
        return False, "N10 position not properly substituted"

    # Check for dimethyl substitution on benzene ring
    # Allow for some flexibility in the exact position
    dimethyl_pattern = Chem.MolFromSmarts("[CH3]~c1~c([CH3])~c~c2")
    if not mol.HasSubstructMatch(dimethyl_pattern):
        return False, "Missing characteristic methyl groups"

    # Additional check for aromaticity of benzene ring
    aromatic_pattern = Chem.MolFromSmarts("c1cccc2")
    if not mol.HasSubstructMatch(aromatic_pattern):
        return False, "Benzene ring not aromatic"

    # Verify overall connectivity and oxidation state
    # Count nitrogens (should have at least 4)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 4:
        return False, "Insufficient number of nitrogen atoms"

    return True, "Contains dimethylisoalloxazine core with N10 substitution"

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