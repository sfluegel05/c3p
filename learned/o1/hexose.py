"""
Classifies: CHEBI:18133 hexose
"""
"""
Classifies: CHEBI:18133 hexose
"""

from rdkit import Chem

def is_hexose(smiles: str):
    """
    Determines if a molecule is a hexose based on its SMILES string.
    A hexose is any six-carbon monosaccharide which in its linear form contains either an aldehyde group at position 1 (aldohexose)
    or a ketone group at position 2 (ketohexose).

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

    # Remove hydrogens for accurate atom counting
    mol = Chem.RemoveHs(mol)

    # Check for invalid atoms (exclude metals and halogens)
    allowed_atomic_nums = {6, 7, 8}  # Allow C, N, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Molecule contains invalid atom type {atom.GetSymbol()}"

    # Define SMARTS patterns for pyranose and furanose rings
    pyranose_pattern = Chem.MolFromSmarts("C1[C@H]([O,C])[C@H]([O,C])[C@@H]([O,C])[C@H]([O,C])O1")
    furanose_pattern = Chem.MolFromSmarts("C1[C@H]([O,C])[C@@H]([O,C])[C@H]([O,C])O1")

    has_pyranose = mol.HasSubstructMatch(pyranose_pattern)
    has_furanose = mol.HasSubstructMatch(furanose_pattern)

    # Check for linear form with aldehyde or ketone
    aldehyde_pattern = Chem.MolFromSmarts("O=CC([O,C])([O,C])[C@H]([O,C])[C@H]([O,C])CO")
    ketone_pattern = Chem.MolFromSmarts("O=C([C@H]([O,C])CO)[C@H]([O,C])[C@H]([O,C])CO")

    has_aldehyde = mol.HasSubstructMatch(aldehyde_pattern)
    has_ketone = mol.HasSubstructMatch(ketone_pattern)

    # Check if molecule contains a monosaccharide unit (no glycosidic bonds)
    glycosidic_bond_pattern = Chem.MolFromSmarts("O[C;R][C;!R]")
    if mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "Molecule appears to have a glycosidic bond, may not be a monosaccharide"

    # Identify hexose core structure
    hexose_pattern = Chem.MolFromSmarts("""
    [$([C@@H]1([O])[C@H]([O])[C@H]([O])[C@@H]([O])[C@H]1[O])],  # Pyranose form
    $([C@@H]1([O])[C@H]([O])[C@@H]([O])[C@H]1[O])               # Furanose form
    ]""")
    has_hexose_core = mol.HasSubstructMatch(hexose_pattern)

    # Final evaluation
    if has_hexose_core or has_pyranose or has_furanose or has_aldehyde or has_ketone:
        return True, "Molecule contains hexose core structure"
    else:
        return False, "Molecule does not match hexose structural criteria"