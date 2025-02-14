"""
Classifies: CHEBI:35366 fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid(smiles: str):
    """
    Determines if a molecule is a fatty acid based on its SMILES string.
    A fatty acid is any aliphatic monocarboxylic acid derived from or contained
    in esterified form in an animal or vegetable fat, oil or wax.
    Natural fatty acids commonly have a chain of 4 to 28 carbons (usually unbranched
    and even-numbered), which may be saturated or unsaturated.
    By extension, the term is sometimes used to embrace all acyclic aliphatic carboxylic acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts('[CX3](=O)[OX1H]')
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    num_carboxylic_acids = len(carboxylic_acid_matches)
    if num_carboxylic_acids == 0:
        return False, "No carboxylic acid group found"
    if num_carboxylic_acids > 1:
        return False, f"Found {num_carboxylic_acids} carboxylic acid groups, expected only one"

    # Check for aromatic rings
    if mol.GetRingInfo().NumAromaticRings() > 0:
        return False, "Contains aromatic rings, not typical for fatty acids"

    # Count the number of carbon atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 4:
        return False, f"Too few carbon atoms ({num_carbons}), not typical for fatty acids"

    # Calculate fraction of sp3 carbons
    num_sp3_carbons = sum(
        1 for atom in mol.GetAtoms()
        if atom.GetAtomicNum() == 6 and atom.GetHybridization() == Chem.HybridizationType.SP3
    )
    fraction_sp3_carbons = num_sp3_carbons / num_carbons if num_carbons > 0 else 0

    # Check if molecule is primarily aliphatic
    if fraction_sp3_carbons < 0.3:
        # Allow unsaturation but ensure aliphatic character
        return False, "Low aliphatic character, molecule may not be a fatty acid"

    return True, "Molecule is a monocarboxylic acid with aliphatic character"