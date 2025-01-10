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
        
    # SMARTS pattern for beta-D-glucose core
    # [OH1,O][C@H]1O[C@H](O)[C@@H](O)[C@H](O)[C@@H]1CO
    # Note: The pattern looks for the beta-D-glucose with correct stereochemistry
    glucose_pattern = Chem.MolFromSmarts('[OX2,OH1][C@H]1O[C@H](O)[C@@H](O)[C@H](O)[C@@H]1CO')
    
    # Alternative pattern that's more flexible with substituents
    glucose_pattern_2 = Chem.MolFromSmarts('[OX2][C@H]1O[C@H]([!H])[C@@H]([OX2])[C@H]([OX2])[C@@H]1CO')
    
    if not (mol.HasSubstructMatch(glucose_pattern) or mol.HasSubstructMatch(glucose_pattern_2)):
        return False, "No beta-D-glucose moiety found with correct stereochemistry"
    
    # Check for glycosidic bond at C1 (anomeric carbon)
    # This is implicit in the match, but we can add additional checks
    
    # Count pyranose rings (6-membered rings containing oxygen)
    ring_pattern = Chem.MolFromSmarts('O1CCCCC1')
    if not mol.HasSubstructMatch(ring_pattern):
        return False, "No pyranose ring found"
    
    # Additional validation: Check for presence of multiple OH groups
    oh_pattern = Chem.MolFromSmarts('O[H]')
    oh_matches = len(mol.GetSubstructMatches(oh_pattern))
    if oh_matches < 3:  # Beta-D-glucose should have multiple OH groups
        return False, "Insufficient hydroxyl groups for beta-D-glucose"
    
    return True, "Contains beta-D-glucose moiety with correct stereochemistry and glycosidic bond"

def test_examples():
    """Test function with some example molecules"""
    examples = [
        "O1[C@@H]([C@H]([C@@H]([C@H]([C@@H]1OC=2C=CC(=CC2)[N+]([O-])=O)O)O)O)CO",  # 4-nitrophenyl beta-D-glucoside
        "OC[C@H]1O[C@@H](Oc2cc3ccc(=O)oc3cc2O)[C@H](O)[C@@H](O)[C@@H]1O",  # esculin
        "CC(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O1"  # Not a beta-D-glucoside
    ]
    
    for smi in examples:
        result, reason = is_beta_D_glucoside(smi)
        print(f"SMILES: {smi}")
        print(f"Is beta-D-glucoside: {result}")
        print(f"Reason: {reason}\n")

if __name__ == "__main__":
    test_examples()