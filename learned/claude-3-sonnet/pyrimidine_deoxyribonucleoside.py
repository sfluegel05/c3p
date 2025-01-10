"""
Classifies: CHEBI:68472 pyrimidine deoxyribonucleoside
"""
"""
Classifies: pyrimidine deoxyribonucleoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_pyrimidine_deoxyribonucleoside(smiles: str):
    """
    Determines if a molecule is a pyrimidine deoxyribonucleoside based on its SMILES string.
    A pyrimidine deoxyribonucleoside consists of a pyrimidine base connected to a 
    deoxyribose sugar via an N-glycosidic bond.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a pyrimidine deoxyribonucleoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for deoxyribose sugar pattern
    # [OH0] means oxygen can't have H attached (ring oxygen)
    # The [CH2] ensures one position lacks an OH (2' deoxy position)
    sugar_pattern = Chem.MolFromSmarts("[OH0;R1]1[CH2;R1][CH1;R1][CH1;R1][CH1;R1]1")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No deoxyribose sugar ring found"

    # Check for pyrimidine base - more general pattern
    # Matches 6-membered ring with two nitrogens in 1,3 positions
    # [nX3] means 3-connected nitrogen (in ring)
    # [#6] means any carbon
    pyrimidine_pattern = Chem.MolFromSmarts("[nX3]1[#6][nX3][#6][#6][#6]1")
    if not mol.HasSubstructMatch(pyrimidine_pattern):
        return False, "No pyrimidine base found"

    # Check for N-glycosidic bond between sugar and base
    # More general pattern that just ensures sugar oxygen and base nitrogen are connected
    # through the right carbons
    glycosidic_pattern = Chem.MolFromSmarts("[OH0;R1]1[CH2][CH1][CH1][CH1]1[NH0;R1]2[#6][#6][#6][#6][#6]2")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No N-glycosidic bond between sugar and base"

    # Verify presence of at least one carbonyl group (characteristic of pyrimidine bases)
    carbonyl_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return False, "Missing characteristic carbonyl group"

    # Check for primary alcohol (5' position)
    primary_oh_pattern = Chem.MolFromSmarts("[CH2]O")
    if not mol.HasSubstructMatch(primary_oh_pattern):
        return False, "Missing 5' hydroxyl group"

    # Check for at least one hydroxyl on sugar (3' position)
    sugar_oh_pattern = Chem.MolFromSmarts("[CH1;R1][OH1,O[C,H]]")
    if not mol.HasSubstructMatch(sugar_oh_pattern):
        return False, "Missing sugar hydroxyl groups"

    # Additional check to ensure the molecule isn't too complex
    # Most pyrimidine nucleosides should have 2 rings (sugar + base)
    rings = mol.GetRingInfo().NumRings()
    if rings > 3:  # Allow up to 3 rings to account for some modifications
        return False, "Structure too complex - too many rings"

    return True, "Contains pyrimidine base connected to deoxyribose via N-glycosidic bond"