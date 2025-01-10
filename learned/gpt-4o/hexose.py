"""
Classifies: CHEBI:18133 hexose
"""
from rdkit import Chem

def is_hexose(smiles: str):
    """
    Determines if a molecule is a hexose based on its SMILES string.
    A hexose is defined as a six-carbon monosaccharide which in its linear form
    contains either an aldehyde group at position 1 (aldohexose) or a ketone group at position 2 (ketohexose),
    or its cyclic form matches common cyclic arrangements like pyranoses or furanoses with ether linkages.

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

    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 6:
        return False, f"Expected 6 carbon atoms, found {c_count}"

    # Define linear aldehyde or ketone patterns at appropriate positions
    aldehyde_pattern = Chem.MolFromSmarts("O=CC([OH1])(C)C")
    ketone_pattern = Chem.MolFromSmarts("O=C(C)C([OH1])([C])C")

    # Define cyclic sugar ring structures (pyranose, furanose) with ether linkages
    pyranose_pattern = Chem.MolFromSmarts("C1(COC1)O")
    furanose_pattern = Chem.MolFromSmarts("C1(CO)CO1")

    # Check for either linear aldehyde/ketone or cyclic forms
    is_aldohexose = mol.HasSubstructMatch(aldehyde_pattern)
    is_ketohexose = mol.HasSubstructMatch(ketone_pattern)
    is_cyclic = mol.HasSubstructMatch(pyranose_pattern) or mol.HasSubstructMatch(furanose_pattern)

    if is_aldohexose:
        return True, "Contains aldohexose group"
    elif is_ketohexose:
        return True, "Contains ketohexose group"
    elif is_cyclic:
        return True, "Contains cyclic hexose structure"
    else:
        return False, "Does not contain aldehyde or ketone group in the right position or cyclic hexose structure"