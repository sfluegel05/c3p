"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
"""
Classifies: beta-D-glucoside compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_D_glucoside(smiles: str):
    """
    Determines if a molecule contains a beta-D-glucoside moiety.
    Beta-D-glucoside has specific stereochemistry:
    - Beta configuration at anomeric carbon (C1)
    - D-configuration for the sugar
    - Pyranose form (6-membered ring)
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a beta-D-glucoside, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for beta-D-glucopyranoside core
    # [O!H] at C1 represents the glycosidic bond
    # Specific stereochemistry at each carbon
    glucose_pattern = Chem.MolFromSmarts('[O!H][C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O1')
    
    # Alternative pattern with more flexible substitution
    glucose_pattern_2 = Chem.MolFromSmarts('[O!H][C@@H]1[C@H]([OH0,OH1])[C@@H]([OH0,OH1])[C@H]([OH0,OH1])[C@@H](CO)O1')
    
    # Pattern specifically for beta configuration at anomeric carbon
    beta_pattern = Chem.MolFromSmarts('[O!H][C@@H]1O[C@H]')
    
    if not (mol.HasSubstructMatch(glucose_pattern) or mol.HasSubstructMatch(glucose_pattern_2)):
        return False, "No beta-D-glucose moiety found with correct stereochemistry"
    
    if not mol.HasSubstructMatch(beta_pattern):
        return False, "Anomeric carbon does not have beta configuration"
    
    # Check for pyranose ring
    ring_pattern = Chem.MolFromSmarts('O1[C][C][C][C][C]1')
    if not mol.HasSubstructMatch(ring_pattern):
        return False, "No pyranose ring found"
    
    # Verify presence of hydroxyl or substituted oxygen groups
    oh_pattern = Chem.MolFromSmarts('[OX2][C]')
    oh_matches = len(mol.GetSubstructMatches(oh_pattern))
    if oh_matches < 4:  # Need at least 4 oxygens (including ring O)
        return False, "Insufficient oxygen-containing groups for beta-D-glucose"
    
    # Additional check for CH2OH group
    ch2oh_pattern = Chem.MolFromSmarts('[CH2]O[H,C,O]')
    if not mol.HasSubstructMatch(ch2oh_pattern):
        return False, "Missing CH2OH group characteristic of glucose"

    # If all checks pass, it's a beta-D-glucoside
    return True, "Contains beta-D-glucose moiety with correct stereochemistry and beta glycosidic bond"

def test_examples():
    """Test function with example molecules"""
    examples = [
        ("CC(CCC[C@@H]1[C@H](CCC1=O)CC(O)=O)O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O",
         "(-)-11-hydroxy-9,10-dihydrojasmonic acid 11-beta-D-glucoside"),
        ("OC[C@H]1O[C@@H](Oc2cc3ccc(=O)oc3cc2O)[C@H](O)[C@@H](O)[C@@H]1O",
         "esculin"),
        ("O1[C@@H]([C@H]([C@@H]([C@H]([C@@H]1OC=2C=CC(=CC2)[N+]([O-])=O)O)O)O)CO",
         "4-nitrophenyl beta-D-glucoside")
    ]
    
    for smi, name in examples:
        result, reason = is_beta_D_glucoside(smi)
        print(f"Name: {name}")
        print(f"SMILES: {smi}")
        print(f"Is beta-D-glucoside: {result}")
        print(f"Reason: {reason}\n")

if __name__ == "__main__":
    test_examples()