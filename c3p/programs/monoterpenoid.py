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

    # Define complex potential patterns indicative of monoterpenoids
    terpene_patterns = [
        Chem.MolFromSmarts("CC(C)C=CC=CCC(C)C"),  # Linear monoterpene example
        Chem.MolFromSmarts("C1CCC(C(C1)C)C"),     # Basic cyclic monoterpene
        Chem.MolFromSmarts("C1C2(C)CCC2CCC1"),    # Bicyclic monoterpenoid
    ]

    # Check for monoterpenoid indicative functional groups, common terpenoid oxygens
    functional_group_patterns = [
        Chem.MolFromSmarts("[CX4][OH]"),  # Alcohol groups
        Chem.MolFromSmarts("[CX3](=O)[OH]"),  # Carboxylic acid groups
        Chem.MolFromSmarts("[CX3](=O)[O][CX4]"),  # Ester groups
    ]

    # Carbon count expected range (allow moderate range around C10 due to modifications)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Check patterns and count
    if ((any(mol.HasSubstructMatch(pattern) for pattern in terpene_patterns) and 8 <= c_count <= 12) or
        (any(mol.HasSubstructMatch(fg_pattern) for fg_pattern in functional_group_patterns) and 8 <= c_count <= 12)):
        return True, "Structure matches typical patterns seen in monoterpenoids with functional groups"

    return False, "No typical monoterpenoid structural features detected"