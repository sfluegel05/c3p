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
    flavanone_core = Chem.MolFromSmarts("[#6]1-[#6]-[#6](=[O])-c2c(-[#6]1)cccc2")
    if flavanone_core is None:
        return False, "Invalid SMARTS pattern for flavanone core"
    
    if not mol.HasSubstructMatch(flavanone_core):
        return False, "No flavanone core structure found"

    # Check for saturated C2-C3 bond (to distinguish from flavones)
    # Match the specific C2-C3 single bond pattern
    c2c3_bond = Chem.MolFromSmarts("O=C1C-C-c2ccccc2O1")
    if c2c3_bond is None:
        return False, "Invalid SMARTS pattern for C2-C3 bond"
        
    if not mol.HasSubstructMatch(c2c3_bond):
        return False, "C2-C3 bond must be saturated (not a flavone)"

    # Check for 4'-hydroxy group
    # This pattern specifically looks for the B-ring with 4'-OH
    # The connection point to the C-ring is specified
    b_ring_4_oh = Chem.MolFromSmarts("[#6]-[#6]1:[#6]:[#6]:[#6](O):[#6]:[#6]:1")
    if b_ring_4_oh is None:
        return False, "Invalid SMARTS pattern for 4'-hydroxy group"
        
    if not mol.HasSubstructMatch(b_ring_4_oh):
        return False, "Missing hydroxy group at 4' position"

    # Additional check to ensure the 4'-OH is properly positioned relative to the core
    complete_pattern = Chem.MolFromSmarts("O=C1CC(c2ccc(O)cc2)Oc2ccccc12")
    if complete_pattern is None:
        return False, "Invalid SMARTS pattern for complete structure"
        
    if not mol.HasSubstructMatch(complete_pattern):
        return False, "4'-hydroxy group not properly positioned relative to core"

    return True, "Molecule is a 4'-hydroxyflavanone"

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