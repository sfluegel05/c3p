"""
Classifies: CHEBI:25409 monoterpenoid
"""
from rdkit import Chem

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.
    A monoterpenoid generally inherits features from a C10 monoterpene backbone, potentially rearranged or modified.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a monoterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Relax the carbon count criteria
    if not (5 <= c_count <= 20):
        return False, f"Carbon count {c_count} not typical for some monoterpenoids but including possible highly modified structures"

    # Expanded set of structural patterns
    structural_patterns = [
        "C1=CC=CC=C1",  # Benzene - to accommodate aromatic transformations
        "C1CCC(CC1)C",  # Cyclohexane - common in monoterpenoids
        "C1CCC(CC1)C=O", # Cyclic ketone patterns
        "C1C=C(C)CC1",  # Cyclic structures often found in monoterpene derivatives
    ]

    for pattern in structural_patterns:
        pattern_mol = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(pattern_mol):
            return True, "Contains structural features typical of modified monoterpene skeletons"

    # Check for presence of common monoterpenoid functional groups, extending previous attempts
    functional_group_patterns = [
        "[CX3](=O)[O,C]",  # Ketone and possible ester
        "[OX2H]",          # Alcohols
        "[NX3;H2,H1;!$(NC=O)]",  # Amines avoiding amides, commonly found in more complex rearrangements
    ]
    
    for fg in functional_group_patterns:
        fg_mol = Chem.MolFromSmarts(fg)
        if mol.HasSubstructMatch(fg_mol):
            return True, "Contains functional groups common to monoterpenoids"

    return False, "Does not have identifiable monoterpenoid features"