"""
Classifies: CHEBI:26935 tetraterpenoid
"""
from rdkit import Chem

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.
    A tetraterpenoid typically has ~40 carbon atoms with conjugated double bonds,
    and may contain cyclic structures or certain functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a tetraterpenoid, False otherwise
        str: Reason for classification or rejection
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check approximate carbon count, allowing some flexibility
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 35 or c_count > 45:
        return False, f"Unsupported carbon count: {c_count} (expected ~40)"

    # Check for sufficient conjugated double bond systems
    polyene_pattern = Chem.MolFromSmarts("C=CC=C")
    polyene_matches = mol.GetSubstructMatches(polyene_pattern)
    if len(polyene_matches) < 3:
        return False, "Insufficient conjugated double-bond systems"

    # Check for typical functional groups like ketones and hydroxyls often found in carotenoids
    functional_group_detected = False
    ketone_pattern = Chem.MolFromSmarts("C=O")
    if mol.HasSubstructMatch(ketone_pattern):
        functional_group_detected = True

    alcohol_pattern = Chem.MolFromSmarts("[OH]")
    if mol.HasSubstructMatch(alcohol_pattern):
        functional_group_detected = True

    if functional_group_detected:
        return True, "Contains functional groups typical for carotenoids (a subset of tetraterpenoids)"

    # Final rule: consider as tetraterpenoid if it meets multiple criteria
    return True, "Meets criteria for a tetraterpenoid, significant polyene chain, and typical functional attribute"

# Example usage
example_smiles = "O=C1C(=C(/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/[C@@H](OC)C(O)(C)C)\\C)\\C)/C)/C)C(C)(C)CC1)C"
is_tetraterpenoid(example_smiles)