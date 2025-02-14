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
    if not (30 <= c_count <= 50):
        return False, f"Expected around 40 carbon atoms, found {c_count}"

    # Check for a more generalized extended conjugated system
    polyene_patterns = [
        Chem.MolFromSmarts("C=CC=CC=CC=C"),  # Existing pattern
        Chem.MolFromSmarts("C=CC=CC=CC=CC=C"),  # Longer conjugation
        Chem.MolFromSmarts("C=C(C=C)C=C")  # Alternate conjugation form
    ]
    found_polyene = any(mol.HasSubstructMatch(pat) for pat in polyene_patterns)
    if not found_polyene:
        return False, "No extended conjugated polyene system found"

    # Check for presence of a diverse set of functional groups
    functional_patterns = [
        Chem.MolFromSmarts("[OH]"),  # Hydroxyl group
        Chem.MolFromSmarts("[C]=O"),  # Carbonyl group
        Chem.MolFromSmarts("[OX2]"),  # Generic oxygen presence, covers hydroxyl and epoxide
        Chem.MolFromSmarts("O"),  # Any oxygen atom
    ]
    found_functional_groups = any(mol.HasSubstructMatch(pat) for pat in functional_patterns)

    # Accept if either polyene system or functional groups are found satisfactorily
    if not found_functional_groups and not found_polyene:
        return False, "Neither functional group nor extended conjugation system found"
    
    return True, "Structure consistent with a tetraterpenoid: extended polyene and/or functional group present"