"""
Classifies: CHEBI:25409 monoterpenoid
"""
from rdkit import Chem

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.
    A monoterpenoid typically has a carbon skeleton derived from a monoterpene
    with around 10 carbons, potentially rearranged or modified.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a monoterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Monoterpenoids generally have around 10 carbons, with flexibility in structural arrangement
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 8 or c_count > 15:
        return False, f"Carbon count {c_count} not typical for monoterpenoid, despite potential rearrangement"

    # Recognize structural flexibility of isoprene-based backbones:
    # e.g., consider various isoprenoid arrangements
    isoprene_like_patterns = [
        Chem.MolFromSmarts("C(C)(C)C(C)=C"),  # Head-to-tail arrangements
        Chem.MolFromSmarts("C=C(C)C(C)"),  # Tail-to-tail arrangements
        Chem.MolFromSmarts("C(C)C(C)=C(C)"),  # C10 branched structures
    ]
    
    if not any(mol.HasSubstructMatch(pattern) for pattern in isoprene_like_patterns):
        return False, "No sufficient isoprene-derived backbone found, not typical for monoterpenoid"

    # Check for presence of common functional groups in monoterpenoids
    common_functional_patterns = {
        "Alcohol": Chem.MolFromSmarts("[OX2H]"),
        "Ketone": Chem.MolFromSmarts("C(=O)[C]"),
        "Ester": Chem.MolFromSmarts("COC(=O)"),
        "Oxide": Chem.MolFromSmarts("COC")
    }

    found_functionalities = [name for name, pattern in common_functional_patterns.items() if mol.HasSubstructMatch(pattern)]
    if found_functionalities:
        return True, f"Contains isoprene-like backbone and functional group(s) {', '.join(found_functionalities)}, consistent with monoterpenoid"

    return True, "Contains isoprene-like backbone without common terminal functional groups but consistent with monoterpenoid structure"