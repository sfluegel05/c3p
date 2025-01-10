"""
Classifies: CHEBI:18133 hexose
"""
from rdkit import Chem

def is_hexose(smiles: str):
    """
    Determines if a molecule is a hexose based on its SMILES string.
    Hexoses are six-carbon monosaccharides which may exist as aldohexoses or ketohexoses 
    in linear forms, and have common cyclic forms as pyranoses or furanoses.
  
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
    
    # Check for exactly 6 carbons
    carbon_pattern = Chem.MolFromSmarts("[#6]")
    c_count = len(mol.GetSubstructMatches(carbon_pattern))
    if c_count != 6:
        return False, f"Incorrect number of carbons for hexose: found {c_count}, expected exactly 6"
    
    # Aldohexose patterns: check for aldehyde group at end of chain
    aldohexose_linear_pattern = Chem.MolFromSmarts("[CX3H1](=O)[C][O;H1]")
    if mol.HasSubstructMatch(aldohexose_linear_pattern):
        return True, "Contains aldehyde group indicating aldohexose"
    
    # Ketohexose patterns: check for ketone group within chain
    ketohexose_linear_pattern = Chem.MolFromSmarts("[CX4][CX3](=O)[CX4][O;H1]")
    if mol.HasSubstructMatch(ketohexose_linear_pattern):
        return True, "Contains ketone group indicating ketohexose"

    # Cyclic aldohexoses and ketohexoses (pyranoses and furanoses)
    pyranose_pattern = Chem.MolFromSmarts("C1(CO)OC[C;R2][C;R2]OC1")
    furanose_pattern = Chem.MolFromSmarts("C1(CO)O[C;R2][C;R2]OC1")
    
    if mol.HasSubstructMatch(pyranose_pattern) or mol.HasSubstructMatch(furanose_pattern):
        return True, "Contains cyclic form indicating hexose"

    return False, "Does not match hexose structure criteria"