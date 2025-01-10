"""
Classifies: CHEBI:24128 furanocoumarin
"""
"""
Classifies: CHEBI:47835 furanocoumarin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_furanocoumarin(smiles: str):
    """
    Determines if a molecule is a furanocoumarin based on its SMILES string.
    A furanocoumarin consists of a furan ring fused with a coumarin core.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a furanocoumarin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More flexible coumarin core patterns that account for substitutions
    coumarin_patterns = [
        Chem.MolFromSmarts("c1cc2c(cc1)OC(=O)CC2"), # Basic coumarin
        Chem.MolFromSmarts("c1cc2c(cc1)OC(=O)C=C2"), # Basic coumarin with double bond
        Chem.MolFromSmarts("[c,C]1[c,C][c,C]2[c,C]([c,C][c,C]1)OC(=O)[C,c][C,c]2"), # More general pattern
    ]
    
    has_coumarin = any(mol.HasSubstructMatch(pat) for pat in coumarin_patterns)
    if not has_coumarin:
        return False, "No coumarin core found"

    # More flexible furan patterns
    furan_patterns = [
        Chem.MolFromSmarts("c1cocc1"), # Basic furan
        Chem.MolFromSmarts("C1COC=C1"), # Dihydrofuran
        Chem.MolFromSmarts("[c,C]1[c,C]o[c,C][c,C]1"), # More general pattern
    ]
    
    has_furan = any(mol.HasSubstructMatch(pat) for pat in furan_patterns)
    if not has_furan:
        return False, "No furan ring found"

    # Patterns for different types of fusion
    fusion_patterns = [
        # Linear furanocoumarins (including variations)
        Chem.MolFromSmarts("c1c2c(ccc1)OC(=O)C=C2c1cocc1"),
        Chem.MolFromSmarts("c1c2c(ccc1)OC(=O)CC2C1COC=C1"),
        Chem.MolFromSmarts("c1c2c(ccc1)OC(=O)C=C2C1=COC=C1"),
        # Angular furanocoumarins (including variations)
        Chem.MolFromSmarts("c1c2cocc2cc2oc(=O)ccc12"),
        Chem.MolFromSmarts("C1c2cocc2cc2oc(=O)ccc12"),
        # More general fusion patterns
        Chem.MolFromSmarts("[c,C]1[c,C]2[c,C]o[c,C][c,C]2[c,C][c,C]2oc(=O)[c,C][c,C][c,C]12"),
        Chem.MolFromSmarts("[c,C]1[c,C]([c,C]2[c,C][c,C]c(=O)o[c,C]2)[c,C]2o[c,C][c,C]2[c,C]1")
    ]

    has_proper_fusion = any(mol.HasSubstructMatch(pat) for pat in fusion_patterns)
    if not has_proper_fusion:
        return False, "Furan ring not properly fused to coumarin core"

    # Ring count check (should have at least 3 rings)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 3:
        return False, "Insufficient number of rings"

    # Determine fusion type
    fusion_type = []
    if any(mol.HasSubstructMatch(Chem.MolFromSmarts(p)) for p in [
        "c1c2c(ccc1)OC(=O)C=C2c1cocc1",
        "c1c2c(ccc1)OC(=O)CC2C1COC=C1"
    ]):
        fusion_type.append("linear")
    if any(mol.HasSubstructMatch(Chem.MolFromSmarts(p)) for p in [
        "c1c2cocc2cc2oc(=O)ccc12",
        "C1c2cocc2cc2oc(=O)ccc12"
    ]):
        fusion_type.append("angular")

    fusion_desc = f"({', '.join(fusion_type)} fusion)" if fusion_type else "(fusion type undetermined)"
    return True, f"Contains furanocoumarin core {fusion_desc}"

def test_furanocoumarin():
    """Test function with known furanocoumarins"""
    test_cases = [
        ("O1C(CC=2C1=C3C(OC(=O)C=C3CCC)=C(C2O)C(=O)C(CC)C)C(O)(C)C", True),  # Cycloneomammein
        ("CC(C)(O)C1Cc2cc3cc(c(=O)oc3cc2O1)C(C)(C)C=C", True),  # Heliettin
        ("Cc1cc2ccc(=O)oc2cc1O", False),  # Simple coumarin
        ("c1ccoc1", False),  # Just furan
    ]
    
    for smiles, expected in test_cases:
        result, reason = is_furanocoumarin(smiles)
        print(f"SMILES: {smiles}")
        print(f"Expected: {expected}, Got: {result}")
        print(f"Reason: {reason}\n")

if __name__ == "__main__":
    test_furanocoumarin()