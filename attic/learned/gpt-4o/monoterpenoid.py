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

    # Monoterpenoids typically have around 10 carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 8 or c_count > 20:  # Expand carbon count range slightly
        return False, f"Carbon count {c_count} not typical for monoterpenoid"

    # Recognize various isoprene-based backbones
    # Consider rearrangements and cyclizations found in monoterpenoids
    isoprene_like_patterns = [
        Chem.MolFromSmarts("C(C)(C)C(C)=C"),  # Head-to-tail
        Chem.MolFromSmarts("C=C(C)C(C)"),  # Tail-to-tail
        Chem.MolFromSmarts("C(C)C(C)=C(C)"),  # C10 branched
        Chem.MolFromSmarts("C1=CC=CC=C1")  # Aromatic rings
    ]
    
    # Check for isoprene-like substructures
    if not any(mol.HasSubstructMatch(pattern) for pattern in isoprene_like_patterns):
        return False, "No sufficient isoprene-derived backbone found, not typical for monoterpenoid"

    # Identify presence of common functional groups
    common_functional_patterns = {
        "Alcohol": Chem.MolFromSmarts("[OX2H]"),
        "Ketone": Chem.MolFromSmarts("C(=O)[C]"),
        "Ester": Chem.MolFromSmarts("COC(=O)"),
        "Oxide": Chem.MolFromSmarts("COC"),
        "Aldehyde": Chem.MolFromSmarts("C=O"),
        "Ether": Chem.MolFromSmarts("C-O-C")
    }

    found_functionalities = [name for name, pattern in common_functional_patterns.items() if mol.HasSubstructMatch(pattern)]
    if found_functionalities:
        return True, f"Contains isoprene-like backbone and functional group(s) {', '.join(found_functionalities)}, consistent with monoterpenoid"

    return True, "Contains isoprene-like backbone without common terminal functional groups but consistent with monoterpenoid structure"