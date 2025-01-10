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
    # More general pattern that allows for substitutions
    flavanone_core = Chem.MolFromSmarts("O=C1CC(c2ccccc2)Oc2ccccc12")
    if not mol.HasSubstructMatch(flavanone_core):
        return False, "No flavanone core structure found"

    # Look for B-ring with 4'-hydroxy
    # This pattern matches the B-ring connected to the C2 position with a hydroxy at 4' position
    b_ring_pattern = Chem.MolFromSmarts("[OH]-[c]1[cH,c][cH,c][c]([CH1]-2-[CH2]-C(=O)[c]3[cH,c][cH,c][cH,c][c]3O2)[cH,c][cH,c]1")
    if not mol.HasSubstructMatch(b_ring_pattern):
        return False, "Missing required 4'-hydroxy group"

    # Verify C2-C3 saturation (to distinguish from flavones)
    # The pattern looks for the specific single bond between C2-C3
    c2c3_pattern = Chem.MolFromSmarts("O=C1[CH2][CH1]([c2ccccc2])Oc2ccccc12")
    if not mol.HasSubstructMatch(c2c3_pattern):
        return False, "Not a flavanone (C2-C3 bond must be saturated)"

    # Additional check for ketone group at C4
    ketone_pattern = Chem.MolFromSmarts("O=C1[CH2][CH1]([c2ccccc2])Oc2ccccc12")
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "Missing required ketone group at C4"

    # Verify the presence of both the flavanone core and 4'-OH
    # Count the number of matches to ensure we have the correct structure
    matches = mol.GetSubstructMatches(b_ring_pattern)
    if len(matches) < 1:
        return False, "Structure does not match 4'-hydroxyflavanone pattern"

    return True, "Molecule is a 4'-hydroxyflavanone with correct core structure and 4'-OH group"

def test_examples():
    """Test function with example molecules"""
    examples = [
        # Known 4'-hydroxyflavanones
        ("O[C@@H]1[C@H](Oc2cc(O)cc(O)c2C1=O)c1ccc(O)c(O)c1", True),  # taxifolin
        ("Oc1ccc(cc1)C1CC(=O)c2c(O)cc(O)cc2O1", True),  # naringenin
        ("COc1cc(O)c(C[C@@H](CC=C(C)C)C(C)=C)c2O[C@@H](CC(=O)c12)c1ccc(O)cc1O", True),  # kurarinone
        
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