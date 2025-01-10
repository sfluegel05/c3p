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

    # Look for coumarin core (benzopyrone)
    # Pattern matches benzene fused with α-pyrone
    coumarin_pattern = Chem.MolFromSmarts("c1cc2c(cc1)OC(=O)C=C2")
    if not mol.HasSubstructMatch(coumarin_pattern):
        return False, "No coumarin core found"

    # Look for furan ring pattern
    furan_pattern = Chem.MolFromSmarts("c1cocc1")
    if not mol.HasSubstructMatch(furan_pattern):
        return False, "No furan ring found"

    # Check for linear furanocoumarin pattern (furo[3,2-g]coumarin)
    linear_pattern = Chem.MolFromSmarts("c1cc2c(cc1)OC(=O)C=C2c1cocc1")
    
    # Check for angular furanocoumarin pattern (furo[2,3-h]coumarin)
    angular_pattern = Chem.MolFromSmarts("c1c2c(cc3oc(=O)ccc3c1)occ2")
    
    if not (mol.HasSubstructMatch(linear_pattern) or mol.HasSubstructMatch(angular_pattern)):
        return False, "Furan ring not properly fused to coumarin core"

    # Count rings to ensure we don't have extra fused rings beyond the basic furanocoumarin
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 3:  # Should have at least 3 rings (benzene, pyrone, furan)
        return False, "Insufficient number of rings"
        
    # Basic checks passed - this appears to be a furanocoumarin
    fusion_type = []
    if mol.HasSubstructMatch(linear_pattern):
        fusion_type.append("linear")
    if mol.HasSubstructMatch(angular_pattern):
        fusion_type.append("angular")
        
    return True, f"Contains furanocoumarin core ({', '.join(fusion_type)} fusion)"

def test_furanocoumarin():
    """Test function with known furanocoumarins"""
    test_cases = [
        ("CC(C)=CCOc1c2ccoc2cc2oc(=O)ccc12", True),  # Simple linear furanocoumarin
        ("Cc1cc2ccc(=O)oc2cc1O", False),  # Coumarin but no furan
        ("c1ccoc1", False),  # Just furan
        ("O1C(C(OC)(C)C)CC2=C1C=CC3=C2OC(=O)C=C3", True),  # Angular furanocoumarin
    ]
    
    for smiles, expected in test_cases:
        result, reason = is_furanocoumarin(smiles)
        print(f"SMILES: {smiles}")
        print(f"Expected: {expected}, Got: {result}")
        print(f"Reason: {reason}\n")

if __name__ == "__main__":
    test_furanocoumarin()