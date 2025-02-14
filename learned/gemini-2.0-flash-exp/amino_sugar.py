"""
Classifies: CHEBI:28963 amino sugar
"""
"""
Classifies: Amino sugar
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar based on its SMILES string.
    An amino sugar is a sugar with one or more hydroxyl groups replaced by substituted or unsubstituted amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino sugar, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS for 5 and 6 membered rings with one oxygen and carbons
    sugar_ring_patterns = [
        Chem.MolFromSmarts("[CX4,CX3]1[OX2][CX4,CX3][CX4,CX3][CX4,CX3]1"), #5 membered ring
        Chem.MolFromSmarts("[CX4,CX3]1[OX2][CX4,CX3][CX4,CX3][CX4,CX3][CX4,CX3]1")  # 6 membered ring
    ]

    has_sugar_ring = False
    for pattern in sugar_ring_patterns:
        if mol.HasSubstructMatch(pattern):
            has_sugar_ring = True
            break

    if not has_sugar_ring:
        return False, "No sugar ring structure found"


    # Check for amino groups (NH2, NHR, NR2) directly bonded to or via a linker to the ring
    amino_group_pattern_direct = Chem.MolFromSmarts("[NX3;!H0]-[CX4,CX3]1[OX2][CX4,CX3][CX4,CX3][CX4,CX3]([CX4,CX3])1") #amino directly bonded to a C in a ring
    substituted_amino_pattern_direct = Chem.MolFromSmarts("[NX3;!H0][CX3](=[OX1])-[CX4,CX3]1[OX2][CX4,CX3][CX4,CX3][CX4,CX3]([CX4,CX3])1") #acetamide directly connected
    amino_group_pattern_indirect = Chem.MolFromSmarts("[NX3;!H0]~[CX4,CX3]~[CX4,CX3]1[OX2][CX4,CX3][CX4,CX3][CX4,CX3]([CX4,CX3])1") #amino group via linker to ring

    amino_matches1 = mol.GetSubstructMatches(amino_group_pattern_direct)
    amino_matches2 = mol.GetSubstructMatches(substituted_amino_pattern_direct)
    amino_matches3 = mol.GetSubstructMatches(amino_group_pattern_indirect)


    if not amino_matches1 and not amino_matches2 and not amino_matches3:
       return False, "No amino group substitution found"

    # If both conditions are met, it is an amino sugar
    return True, "Contains a sugar ring with one or more amino group substitutions."