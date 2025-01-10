"""
Classifies: CHEBI:18133 hexose
"""
from rdkit import Chem

def is_hexose(smiles: str):
    """
    Determines if a molecule is a hexose based on its SMILES string.
    A hexose is a six-carbon monosaccharide which in its linear form
    contains an aldehyde group at position 1 (aldohexose) or a ketone group at position 2 (ketohexose),
    or in its cyclic form matches pyranose or furanose structures.
    
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

    # Count carbon atoms to ensure it's a hexose
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 6:
        return False, f"Expected 6 carbon atoms, found {c_count}"

    # Define broader SMARTS patterns for aldohexose and ketohexose
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[CH2][CX4][CX4][CX4][CH2][OH]")  # General aldose pattern
    ketone_pattern = Chem.MolFromSmarts("[CH2][CX3](=O)[CH][CX4][CH2][CH2][OH]")  # General ketose pattern

    # Define more flexible cyclic pyranose and furanose patterns
    pyranose_pattern = Chem.MolFromSmarts("C1OC(O)C(O)C(O)C1")  # Broader pyranose pattern
    furanose_pattern = Chem.MolFromSmarts("C1OC(O)C(O)C1")  # Broader furanose pattern

    # Check for linear aldohexose or ketohexose
    is_aldohexose = mol.HasSubstructMatch(aldehyde_pattern)
    is_ketohexose = mol.HasSubstructMatch(ketone_pattern)

    # Check for cyclic pyranose or furanose structures
    is_pyranose = mol.HasSubstructMatch(pyranose_pattern)
    is_furanose = mol.HasSubstructMatch(furanose_pattern)

    if is_aldohexose:
        return True, "Contains aldohexose structure"
    elif is_ketohexose:
        return True, "Contains ketohexose structure"
    elif is_pyranose or is_furanose:
        return True, "Contains cyclic hexose structure"
    else:
        return False, "Does not match hexose structure in either linear or cyclic forms"