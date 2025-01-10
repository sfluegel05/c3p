"""
Classifies: CHEBI:18133 hexose
"""
from rdkit import Chem

def is_hexose(smiles: str):
    """
    Determines if a molecule is a hexose based on its SMILES string.
    Hexoses can be six-carbon monosaccharides which contain either an aldehyde group or a ketone group in their linear form,
    or can be present in typical cyclic forms as pyranoses or furanoses.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hexose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for exactly 6 carbons in the main framework (some flexibility for cyclic forms)
    carbon_pattern = Chem.MolFromSmarts("[#6]")
    c_count = len(mol.GetSubstructMatches(carbon_pattern))
    if c_count < 6:
        return False, f"Insufficient number of carbons for hexose: found {c_count}, expected at least 6"
    
    # Recognize typical hexose forms:
    # 1. Linear aldohexoses patterns like CHO-CH2-(OH)4
    # 2. Linear ketohexoses patterns like CH2OH-CO-CH-(OH)3
    # 3. Cyclic forms including pyranoses and furanoses (rings with any 5 or 6 atoms containing oxygen)
    
    # Aldohexose patterns
    aldohexose_linear_pattern = Chem.MolFromSmarts("[CX3H](=O)[CX4][OX2H]")
    if mol.HasSubstructMatch(aldohexose_linear_pattern):
        return True, "Contains aldehyde group indicating aldohexose"
    
    # Ketohexose patterns
    ketohexose_linear_pattern = Chem.MolFromSmarts("[CX4H2][CX3](=O)[CX4][OX2H]")
    if mol.HasSubstructMatch(ketohexose_linear_pattern):
        return True, "Contains ketone group indicating ketohexose"

    # Cyclic hexoses (both aldo and keto forms)
    pyranose_pattern = Chem.MolFromSmarts("[C&R1][O&R][C&R][C&R][C&R][O&R1]") # 6-membered ring with O
    furanose_pattern = Chem.MolFromSmarts("[C&R1][O&R][C&R][C&R][O&R1]") # 5-membered ring with O
    
    if mol.HasSubstructMatch(pyranose_pattern) or mol.HasSubstructMatch(furanose_pattern):
        return True, "Contains cyclic form indicating hexose"

    return False, "Does not match hexose structure criteria"