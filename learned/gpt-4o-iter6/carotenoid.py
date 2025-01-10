"""
Classifies: CHEBI:23044 carotenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_carotenoid(smiles: str):
    """
    Determines if a molecule is a carotenoid based on its SMILES string.
    A carotenoid is a tetraterpenoid (C40) with characteristic structural features.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carotenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carbon count around 40 (reflective of tetraterpenoids)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (38 <= c_count <= 42):
        return False, f"Carbon count {c_count} not consistent with tetraterpenoid structure"

    # Count conjugated double bonds (pattern C=C-C=C-...)
    conjugated_pattern = Chem.MolFromSmarts("C=C")
    conjugated_matches = mol.GetSubstructMatches(conjugated_pattern)
    if len(conjugated_matches) < 8:
        return False, f"Conjugation pattern too short; found {len(conjugated_matches)} double bonds"

    # Check for other common structural modifications - needing some functional groups
    carbonyl_oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetDegree() > 1)
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1)

    if carbonyl_oxygen_count == 0 and hydroxyl_count == 0:
        return False, "No oxidation pattern detected (e.g., hydroxyl or carbonyl groups)"

    # Check if saturated rings can be identified
    cyclic_patterns = Chem.MolFromSmarts("C1CCCC1")
    if not mol.HasSubstructMatch(cyclic_patterns):
        return False, "No notable cyclization (rings) detected"

    return True, "Structure consistent with carotenoid (C40 tetraterpenoid) with conjugation and modifications"