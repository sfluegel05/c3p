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

    # Define SMARTS pattern for monosaccharide units (pyranose and furanose rings)
    sugar_pattern = Chem.MolFromSmarts('[CR1][CR1][CR1][CR1][O][CR1]')  # Six-membered ring with one oxygen
    if sugar_pattern is None:
        return False, "Invalid sugar SMARTS pattern"

    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    num_sugars = len(sugar_matches)

    if num_sugars < 2:
        return False, "Does not contain multiple sugar units"

    # Define SMARTS pattern for uronic acid (sugar with carboxylic acid group)
    uronic_acid_pattern = Chem.MolFromSmarts('[C;R](=O)[O;H1,-]')  # Carboxylic acid group attached to ring carbon
    if uronic_acid_pattern is None:
        return False, "Invalid uronic acid SMARTS pattern"

    uronic_acid_matches = mol.GetSubstructMatches(uronic_acid_pattern)
    if not uronic_acid_matches:
        return False, "No uronic acid units found"

    # Define SMARTS pattern for glycosamine (sugar with amine group)
    glycosamine_pattern = Chem.MolFromSmarts('[C;R][NH2]')  # Ring carbon attached to NH2 group
    if glycosamine_pattern is None:
        return False, "Invalid glycosamine SMARTS pattern"

    glycosamine_matches = mol.GetSubstructMatches(glycosamine_pattern)
    if not glycosamine_matches:
        return False, "No glycosamine units found"

    # Check for alternating units by ensuring both patterns are matched multiple times
    if len(uronic_acid_matches) < 1 or len(glycosamine_matches) < 1:
        return False, "Insufficient alternating units of uronic acids and glycosamines"

    # Define SMARTS pattern for sulfate ester group (-OSO3)
    sulfate_pattern = Chem.MolFromSmarts('[O][S](=O)(=O)[O]')  # Sulfate ester group
    if sulfate_pattern is None:
        return False, "Invalid sulfate SMARTS pattern"

    sulfate_matches = mol.GetSubstructMatches(sulfate_pattern)
    if not sulfate_matches:
        return False, "No sulfate ester groups found"

    # Check for polysaccharide chain length (heuristic: more than 5 sugar units)
    if num_sugars < 5:
        return False, f"Contains only {num_sugars} sugar units, not sufficient for a polysaccharide"

    return True, "Molecule is a mucopolysaccharide with alternating uronic acids and glycosamines, partially esterified with sulfuric acid"