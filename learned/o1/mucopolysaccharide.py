"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_mucopolysaccharide(smiles: str):
    """
    Determines if a molecule is a mucopolysaccharide based on its SMILES string.
    Mucopolysaccharides are polysaccharides composed of alternating units from uronic acids and
    glycosamines, and commonly partially esterified with sulfuric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mucopolysaccharide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for sugar rings (pyranose and furanose)
    sugar_pattern = Chem.MolFromSmarts("""
        [$([C@@H]1O[C@H]([C@@H](O)[C@H](O)[C@H]1O])],
        [$([C@@H]1O[C@@H]([C@H](O)[C@@H](O)[C@H]1O])]
    """)
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if len(sugar_matches) < 2:
        return False, "Does not contain multiple sugar units"

    # Define SMARTS patterns for uronic acid (sugar acids with carboxylic group)
    uronic_acid_pattern = Chem.MolFromSmarts("""
        [C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O
        -[#6](=O)[O-]
    """)
    uronic_acid_matches = mol.GetSubstructMatches(uronic_acid_pattern)
    if not uronic_acid_matches:
        return False, "No uronic acid units found"

    # Define SMARTS pattern for glycosamine (sugar with amine group)
    glycosamine_pattern = Chem.MolFromSmarts("""
        [C@@H]1O[C@H](CN)[C@@H](O)[C@H](O)[C@H]1O
    """)
    glycosamine_matches = mol.GetSubstructMatches(glycosamine_pattern)
    if not glycosamine_matches:
        return False, "No glycosamine units found"

    # Check for alternating units by ensuring both patterns are matched multiple times
    if len(uronic_acid_matches) < 1 or len(glycosamine_matches) < 1:
        return False, "Insufficient alternating units of uronic acids and glycosamines"

    # Define SMARTS pattern for sulfate ester group (-OSO3)
    sulfate_pattern = Chem.MolFromSmarts("[$([OX2]S(=O)(=O)[O-])]")
    sulfate_matches = mol.GetSubstructMatches(sulfate_pattern)

    if not sulfate_matches:
        return False, "No sulfate ester groups found"

    # Check for polysaccharide chain length (heuristic: more than 5 sugar units)
    num_sugars = len(sugar_matches)
    if num_sugars < 5:
        return False, f"Contains only {num_sugars} sugar units, not sufficient for a polysaccharide"

    return True, "Molecule is a mucopolysaccharide with alternating uronic acids and glycosamines, partially esterified with sulfuric acid"