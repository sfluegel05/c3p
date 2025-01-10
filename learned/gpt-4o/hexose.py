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

    # Define more accurate linear SMARTS patterns for aldehyde at position 1 and ketone at position 2
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[CH2][CH]([OH])[CH]([OH])[CH]([OH])[CH2][OH]")  # Aldohexose
    ketone_pattern = Chem.MolFromSmarts("[CH2][CX3](=O)[CH]([OH])[CH]([OH])[CH2][CH2][OH]")  # Ketohexose

    # Define patterns for pyranose and furanose rings
    pyranose_pattern = Chem.MolFromSmarts("C1[C@H]([O])[C@@H]([OH])[C@H]([OH])[C@@H]([OH])O1")  # Pyranose
    furanose_pattern = Chem.MolFromSmarts("C1[C@H]([O])[C@@H]([OH])[C@H]([OH])O1")  # Furanose

    # Check for linear aldohexose or ketohexose
    is_aldohexose = mol.HasSubstructMatch(aldehyde_pattern)
    is_ketohexose = mol.HasSubstructMatch(ketone_pattern)

    # Check for cyclic pyranose or furanose structures
    is_pyranose = mol.HasSubstructMatch(pyranose_pattern)
    is_furanose = mol.HasSubstructMatch(furanose_pattern)

    if is_aldohexose:
        return True, "Contains aldohexose group"
    elif is_ketohexose:
        return True, "Contains ketohexose group"
    elif is_pyranose or is_furanose:
        return True, "Contains cyclic hexose structure"
    else:
        return False, "Does not match hexose structure in either linear or cyclic forms"