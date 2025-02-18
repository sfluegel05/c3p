"""
Classifies: CHEBI:25409 monoterpenoid
"""
from rdkit import Chem

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.
    Monoterpenoids are derived from monoterpenes (C10), potentially with rearrangements
    or modifications by the removal or addition of functional groups.

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

    # Define expanded potential patterns indicative of monoterpenoids
    terpene_patterns = [
        Chem.MolFromSmarts("C1=CCC(C=C1)C(C)C"),     # Common cyclic monoterpenoid (p-menthane derivatives)
        Chem.MolFromSmarts("C12CCC(C1)(C)C2"),       # Bicyclic monoterpenoids
        Chem.MolFromSmarts("C1=CC=C(C(C)C)CC1"),     # Aromatic monoterpenoid variations
    ]

    # Check for common monoterpenoid functional groups within expected span
    functional_group_patterns = [
        Chem.MolFromSmarts("[CX4][OH]"),  # Alcohol groups directly attached to terpenoid structure
        Chem.MolFromSmarts("[CX3](=O)[OH]"),  # Carboxylic acid groups on modified terpenoids
        Chem.MolFromSmarts("[CX3](=O)[O][CX4]"),  # Ester groups common in monoterpenoid derivatives
    ]

    # Carbon count expected to be around 10 (allowing slight flexibility for rearrangements and functional groups)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Check patterns and count
    if ((any(mol.HasSubstructMatch(pattern) for pattern in terpene_patterns) and 9 <= c_count <= 11) or
        (all(mol.HasSubstructMatch(fg_pattern) for fg_pattern in functional_group_patterns) and 9 <= c_count <= 11)):
        return True, "Structure matches typical patterns seen in monoterpenoids with functional groups"

    return False, "No typical monoterpenoid structural features detected"