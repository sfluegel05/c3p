"""
Classifies: CHEBI:26935 tetraterpenoid
"""
from rdkit import Chem

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetraterpenoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (35 <= c_count <= 42):
        return False, f"Expected around 40 carbon atoms, found {c_count}"

    # Check for conjugated double bonds pattern
    # Adjusting SMARTS for longer conjugation typically seen in tetraterpenoids
    polyene_pattern = Chem.MolFromSmarts("C=CC=CC=CC=C")
    if not mol.HasSubstructMatch(polyene_pattern):
        return False, "No extended conjugated polyene system found"

    # Check for presence of functional groups (e.g., -OH, =O)
    functional_patterns = [
        Chem.MolFromSmarts("[OH]"),  # Hydroxyl group
        Chem.MolFromSmarts("[C]=O"),  # Carbonyl group
        Chem.MolFromSmarts("[OX2]"),  # Generic oxygen presence, covers hydroxyl and epoxide
    ]
    found_functional_groups = any(mol.HasSubstructMatch(pat) for pat in functional_patterns)
    if not found_functional_groups:
        return False, "No common functional groups (e.g., hydroxyl, carbonyl, epoxide) found"

    return True, "Structure consistent with a tetraterpenoid: extended polyene system and functional group present"