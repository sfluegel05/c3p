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
    
    # Relax the carbon count criteria to ensure highly modified structures are included
    if not (5 <= c_count <= 30):  # Extended upper range to 30
        return False, f"Carbon count {c_count} not typical for common monoterpenoids but may include highly modified structures"

    # Expanded set of structural patterns common in monoterpenoids
    structural_patterns = [
        "C1=CC=CC=C1",          # Benzene - to accommodate aromatic transformations
        "C1CCC(CC1)C",          # Cyclohexane - common in monoterpenoids
        "C1C=C(C)CC1",          # Small cyclic structures in terpene derivatives
        "C(C)=C",               # Acyclic terpenoid characteristic
        "C1C=C(CC1)C",          # Larger cyclic structures in monoterpene derivatives
    ]

    for pattern in structural_patterns:
        pattern_mol = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(pattern_mol):
            return True, "Contains structural features typical of modified monoterpene skeletons"

    # Check for presence of common monoterpenoid functional groups, further extended
    functional_group_patterns = [
        "[CX3](=O)[O,C]",       # Ketone and possible ester groups
        "[OX2H]",               # Alcohols
        "[NX3;H2,H1;!$(NC=O)]", # Amines avoiding amides
        "[OX2][CX3](=O)[CX3]",  # Ester or lactone configurations
        "[CX3](=[OX1])[OX2H0]C",# Ethers in conjunction with carbonyl groups
    ]
    
    for fg in functional_group_patterns:
        fg_mol = Chem.MolFromSmarts(fg)
        if mol.HasSubstructMatch(fg_mol):
            return True, "Contains functional groups common to monoterpenoids"

    return False, "Does not have identifiable monoterpenoid features"