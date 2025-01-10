"""
Classifies: CHEBI:23044 carotenoid
"""
from rdkit import Chem

def is_carotenoid(smiles: str):
    """
    Determines if a molecule is a carotenoid based on its SMILES string.
    A carotenoid is characterized as a tetraterpenoid (C40) with flexibility in structure, conjugation, and functional derivatives.

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
    if c_count < 30 or c_count > 52:
        return False, f"Carbon count {c_count} is not typical for carotenoids"

    # Check for long conjugated C=C system with flexibility in patterns
    conjugation_patterns = [
        "C=C-C=C-C=C",  # extended conjugation
        "C=C-C=C-C=C-C=C",
    ]
    if not any(mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)) for pattern in conjugation_patterns):
        return False, "Lacks extended conjugated system typical of carotenoids"

    # Check common functional groups in carotenoids (possible oxy groups)
    h_oxygen = Chem.MolFromSmarts("[OX2H]")  # Hydroxyl group
    ketone = Chem.MolFromSmarts("C=O")  # Ketone group
    if not mol.HasSubstructMatch(h_oxygen) and not mol.HasSubstructMatch(ketone):
        return False, "Few oxygens and no typical carotenoid functional groups"

    # Carotenoids can have cyclization, check for flexible ring structures
    ring_found = any(mol.GetRingInfo().IsAtomInRingOfSize(atom.GetIdx(), size) for atom in mol.GetAtoms() for size in range(5, 7))
    if c_count >= 36 and not ring_found:
        return False, "Carotenoids generally contain cyclic derivatives"

    return True, "Structure fits the classification of a carotenoid based on core patterns"