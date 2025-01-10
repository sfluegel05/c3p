"""
Classifies: CHEBI:26935 tetraterpenoid
"""
from rdkit import Chem

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.
    A tetraterpenoid typically has 40 carbon atoms with conjugated double bonds,
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

    # Check carbon count
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 40:
        return False, f"Incorrect number of carbons: {c_count} (expected 40)"

    # Check for conjugated double bonds (pattern: CC=CC or similar)
    conjugated_pattern = Chem.MolFromSmarts("C=C-C=C")
    conjugated_matches = mol.GetSubstructMatches(conjugated_pattern)
    if len(conjugated_matches) < 1:
        return False, "No conjugated double bonds found"
    
    # Optional: Check for typical functional groups
    ketone_pattern = Chem.MolFromSmarts("C=O")
    if mol.HasSubstructMatch(ketone_pattern):
        return True, "Contains ketone group, typical for carotenoids (subset of tetraterpenoids)"

    alcohol_pattern = Chem.MolFromSmarts("O")
    if mol.HasSubstructMatch(alcohol_pattern):
        return True, "Contains hydroxyl group, likely to be a tetraterpenoid"

    return True, "Meets criteria for a tetraterpenoid"

# Example usage
example_smiles = "O=C1C(=C(/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/[C@@H](OC)C(O)(C)C)\\C)\\C)/C)/C)C(C)(C)CC1)C"
is_tetraterpenoid(example_smiles)