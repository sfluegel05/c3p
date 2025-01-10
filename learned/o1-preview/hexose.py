"""
Classifies: CHEBI:18133 hexose
"""
from rdkit import Chem

def is_hexose(smiles: str):
    """
    Determines if a molecule is a hexose based on its SMILES string.
    A hexose is any six-carbon monosaccharide which in its linear form contains either
    an aldehyde group at position 1 (aldohexose) or a ketone group at position 2 (ketohexose).
    This function accounts for both linear and cyclic forms.

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

    # Remove salts and separate components, keep only the largest fragment
    mol = Chem.RemoveHs(mol)
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    mol = max(frags, default=mol, key=lambda m: m.GetNumAtoms())

    # Count number of carbons
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons != 6:
        return False, f"Number of carbon atoms ({num_carbons}) is not 6"

    # Check for monosaccharide structure (no glycosidic bonds)
    # Count the number of glycosidic bonds (C-O-C with one carbon being anomeric)
    glyco_bond_pattern = Chem.MolFromSmarts("[#6]-O-[#6]")
    num_glyco_bonds = len(mol.GetSubstructMatches(glyco_bond_pattern))
    if num_glyco_bonds > 1:
        return False, "Molecule has multiple glycosidic bonds, may not be a monosaccharide"

    # Define SMARTS patterns for aldohexoses and ketohexoses in linear and cyclic forms
    # Aldohexose linear form: O=CH-CH(-OH)-CH(-OH)-CH(-OH)-CH(-OH)-CH2OH
    aldohexose_linear = Chem.MolFromSmarts("O=CO[C@H](O)[C@H](O)[C@H](O)CO")
    # Ketohexose linear form: HO-CH2-C(=O)-CH(-OH)-CH(-OH)-CH(-OH)-CH2OH
    ketohexose_linear = Chem.MolFromSmarts("O=CC(O)[C@H](O)[C@H](O)CO")
    # Aldopyranose (6-membered ring)
    aldohexose_cyclic = Chem.MolFromSmarts("C1[C@H](O)[C@H](O)[C@H](O)[C@H](O)O1")
    # Aldofuranose (5-membered ring)
    aldohexose_furanose = Chem.MolFromSmarts("C1[C@H](O)[C@H](O)[C@H](O)O1")
    # Ketopyranose (6-membered ring)
    ketohexose_cyclic = Chem.MolFromSmarts("C1[C@H](O)[C@H](O)[C@H](O)[C@@H](O)C(=O)O1")
    # Ketofuranose (5-membered ring)
    ketohexose_furanose = Chem.MolFromSmarts("C1[C@H](O)[C@H](O)[C@@H](O)C(=O)O1")

    # Verify that the SMARTS patterns are valid
    patterns = [aldohexose_linear, ketohexose_linear, aldohexose_cyclic,
                aldohexose_furanose, ketohexose_cyclic, ketohexose_furanose]
    pattern_names = ["aldohexose linear", "ketohexose linear", "aldohexose pyranose",
                     "aldohexose furanose", "ketohexose pyranose", "ketohexose furanose"]

    for pattern, name in zip(patterns, pattern_names):
        if pattern is None:
            return False, f"SMARTS pattern for {name} is invalid"

    # Match patterns
    if mol.HasSubstructMatch(aldohexose_linear):
        return True, "Molecule is an aldohexose in linear form"
    elif mol.HasSubstructMatch(ketohexose_linear):
        return True, "Molecule is a ketohexose in linear form"
    elif mol.HasSubstructMatch(aldohexose_cyclic):
        return True, "Molecule is an aldohexose in pyranose form"
    elif mol.HasSubstructMatch(aldohexose_furanose):
        return True, "Molecule is an aldohexose in furanose form"
    elif mol.HasSubstructMatch(ketohexose_cyclic):
        return True, "Molecule is a ketohexose in pyranose form"
    elif mol.HasSubstructMatch(ketohexose_furanose):
        return True, "Molecule is a ketohexose in furanose form"

    # As an alternative, check for general hexose pattern
    # Six carbons, multiple hydroxyl groups, and possibly one carbonyl group
    num_oxy = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if num_oxy >= 6:
        return True, "Molecule has six carbons and multiple hydroxyl groups; likely a hexose"

    return False, "Molecule does not match hexose patterns"