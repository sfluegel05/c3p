"""
Classifies: CHEBI:23044 carotenoid
"""
from rdkit import Chem

def is_carotenoid(smiles: str):
    """
    Determines if a molecule is a carotenoid based on its SMILES string.
    A carotenoid is characterized as a tetraterpenoid (C40) with possible structural variations.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a carotenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Count carbon atoms - carotenoids typically have a C40 backbone with flexibility
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 30 or c_count > 50:
        return False, f"Carbon count {c_count} is not typical for carotenoids"

    # Look for extended conjugation: characteristic long chains of conjugated C=C
    conjugated_pattern = Chem.MolFromSmarts("C=C-C=C")
    if not mol.HasSubstructMatch(conjugated_pattern):
        return False, "Lacks extended conjugated system typical of carotenoids"

    # Check for presence of functional groups, which are common in carotenoid derivatives
    # Includes hydroxyl (OH), ketone (C=O), etc.
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if c_count < 36 and o_count == 0:
        return False, "Few carbons and no oxygen; atypical for this class"

    # Carotenoids may have a ring structure derived from cyclizations (flexible check)
    possible_ring_smarts = ["C1CCCC1", "C1CCCCC1", "C1C=CCC1", "C1C=CCCC1"]
    if not any(mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)) for smarts in possible_ring_smarts):
        return False, "No ring or cyclic derivatives detected; inconsistent with typical carotenoids"

    return True, "Structure fits the classification of a carotenoid based on core patterns"