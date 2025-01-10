"""
Classifies: CHEBI:140331 4'-hydroxyflavanones
"""
"""
Classifies: 4'-hydroxyflavanones
Definition: Any hydroxyflavanone having a hydroxy substituent located at position 4'
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_4__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 4'-hydroxyflavanone based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (is_4_hydroxyflavanone, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic flavanone core (chroman-4-one system)
    # More flexible pattern that allows for substitutions
    flavanone_core = Chem.MolFromSmarts("O=C1CC(Ar)Oc2cccc12")
    if not mol.HasSubstructMatch(flavanone_core):
        return False, "No flavanone core structure found"

    # Look for the B-ring with 4'-hydroxy group
    # This pattern specifically identifies the 4' position relative to the attachment point
    b_ring_pattern = Chem.MolFromSmarts("[OH]c1ccc(C2CC(=O)c3ccccc3O2)cc1")
    
    # Alternative pattern that allows for substitutions at other positions
    b_ring_pattern_2 = Chem.MolFromSmarts("[OH]c1cc([c,n])c(C2CC(=O)c3c([*,H])c([*,H])c([*,H])c([*,H])c3O2)c([c,n])1")
    
    if not (mol.HasSubstructMatch(b_ring_pattern) or mol.HasSubstructMatch(b_ring_pattern_2)):
        return False, "Missing required 4'-hydroxy group"

    # Verify it's a flavanone (not a flavone) by checking C2-C3 saturation
    # The pattern looks for the saturated C2-C3 bond
    c2c3_pattern = Chem.MolFromSmarts("O=C1CC([#6])Oc2ccccc12")
    if not mol.HasSubstructMatch(c2c3_pattern):
        return False, "Not a flavanone (C2-C3 bond must be saturated)"

    # Additional validation to ensure proper structure
    # Count the number of carbonyl groups in the correct position
    ketone_pattern = Chem.MolFromSmarts("O=C1CC([#6])Oc2ccccc12")
    if len(mol.GetSubstructMatches(ketone_pattern)) != 1:
        return False, "Incorrect number of ketone groups"

    # Check for any obvious disqualifying features
    # Like an unsaturated C2-C3 bond characteristic of flavones
    flavone_pattern = Chem.MolFromSmarts("O=C1C=C([#6])Oc2ccccc12")
    if mol.HasSubstructMatch(flavone_pattern):
        return False, "Structure appears to be a flavone rather than a flavanone"

    return True, "Molecule is a 4'-hydroxyflavanone with correct core structure and 4'-OH group"

def test_examples():
    """Test function with example molecules"""
    examples = [
        # Known 4'-hydroxyflavanones
        ("COc1cc(O)c(C[C@@H](CC=C(C)C)C(C)=C)c2O[C@@H](CC(=O)c12)c1ccc(O)cc1O", True),  # (2S)-(-)-kurarinone
        ("O[C@@H]1[C@H](Oc2c(O)c(O)cc(O)c2C1=O)c1ccc(O)c(O)c1", True),  # 2,3-dihydrogossypetin
        ("Oc1ccc(cc1)C1CC(=O)c2c(O)cc(O)cc2O1", True),  # naringenin
        
        # Negative examples
        ("COc1cc(O)cc2OC(CC(=O)c12)c1ccccc1", False),  # missing 4'-OH
        ("O=C1CC(c2ccccc2)Oc2ccccc12", False),  # no hydroxy groups
    ]
    
    for smiles, expected in examples:
        result, reason = is_4__hydroxyflavanones(smiles)
        print(f"SMILES: {smiles}")
        print(f"Expected: {expected}, Got: {result}")
        print(f"Reason: {reason}\n")

if __name__ == "__main__":
    test_examples()