"""
Classifies: CHEBI:33916 aldopentose
"""
from rdkit import Chem

def is_aldopentose(smiles: str):
    """
    Determines if a molecule is an aldopentose based on its SMILES string.

    An aldopentose is a pentose (5-carbon sugar) with a potential aldehyde
    group at one end when in linear form. Often exists in cyclic form (furanose 
    or pyranose) within the context of hemiacetal structures.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an aldopentose, False otherwise.
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if there are exactly 5 carbon atoms (pentose sugar)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 5:
        return False, f"Expected 5 carbon atoms, found {c_count}"
    
    # Look for cyclic furanose (5 ring) or pyranose (6 ring) structures typical for sugars
    furanose_pattern = Chem.MolFromSmarts("C1OC(O)C(O)C1")  # 5-member furanose
    pyranose_pattern = Chem.MolFromSmarts("C1OC(O)C(O)C(O)C1")  # 6-member pyranose
    
    if mol.HasSubstructMatch(furanose_pattern) or mol.HasSubstructMatch(pyranose_pattern):
        return True, "Contains cyclic furanose or pyranose structure typical for aldopentoses"
    
    # Check for the presence of an aldehyde group (R-CHO) with adjacent hydroxyl
    aldehyde_pattern = Chem.MolFromSmarts("[#6][CX3H1](=O)[OH]")  # Aldehyde with next carbon hydroxyl
    if mol.HasSubstructMatch(aldehyde_pattern):
        return True, "Contains an appropriate 5-carbon chain with an aldehyde group at one end and hydroxyl groups"

    return False, "Does not match aldopentose structure requirements"