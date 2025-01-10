"""
Classifies: CHEBI:74716 withanolide
"""
from rdkit import Chem

def is_withanolide(smiles: str):
    """
    Determines if a molecule is a withanolide based on its SMILES string.
    A withanolide should typically include a modified steroid skeleton,
    including additional oxygen functionalities and a lactone group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a withanolide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the steroid core with a more flexible SMARTS pattern
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C3CC(C(C(C3)C2)C4CCC5C4(CCC5)O)C1")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid-like structure found"
    
    # Look for lactone groups attached to the steroid D-ring
    lactone_pattern = Chem.MolFromSmarts("C1=CC(=O)OCC1")  # Lactone structure within a ring
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone ring attached to steroid structure found"
        
    # Detect common oxygen functionalities
    oxy_func_patterns = [
        Chem.MolFromSmarts("O"),  # General oxygen presence (hydroxyl, ketone, etc.)
        Chem.MolFromSmarts("C=O"), # Carbonyl group
        Chem.MolFromSmarts("OC"),
        Chem.MolFromSmarts("CO"),
    ]

    # Check that there are several oxygens in various functionalities
    oxy_count = sum(len(mol.GetSubstructMatches(patt)) for patt in oxy_func_patterns)
    if oxy_count < 4:  # Arbitrary threshold for withanolides which are heavily oxygenated
        return False, "Insufficient oxygen functionalities detected"

    return True, "Matches characteristics of withanolide: steroid core, lactone, and oxygen functionalities"