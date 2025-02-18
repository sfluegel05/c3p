"""
Classifies: CHEBI:25409 monoterpenoid
"""
from rdkit import Chem

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.
    Monoterpenoids are derivatives of monoterpenes, containing a 10-carbon skeleton with potential rearrangements.

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

    # Monoterpenoid is characterized by having a 10-carbon backbone
    # including cyclic or acyclic structures derived from isoprene units
    monoterpenoid_patterns = [
        Chem.MolFromSmarts("C1CCCCC1"),     # Monocyclic terpenes (like menthane derivatives)
        Chem.MolFromSmarts("C1CCC(C)CC1"),  # Cyclohexane modified by isoprene units
        Chem.MolFromSmarts("C1=CC=CC=C1"),  # Aromatic rings as part of rearrangement
        Chem.MolFromSmarts("C1CCC(C)C1"),   # Simple cyclopentane
        Chem.MolFromSmarts("C2(CC=CC1C2)C1"),  # Bicyclic structures like pinanes
        Chem.MolFromSmarts("C=C(C)C"),      # Acyclic monoterpenes
    ]

    # Functional groups frequently found in monoterpenoids
    functional_group_patterns = [
        Chem.MolFromSmarts("[OX2H]"),       # Hydroxyl group
        Chem.MolFromSmarts("[CX3]=[OX1]"),  # Carbonyl group
        Chem.MolFromSmarts("[CX3](=O)[OX2H1]"),  # Carboxylic acids
        Chem.MolFromSmarts("O=C(O)C"),      # Esters
        Chem.MolFromSmarts("C=CC"),         # Double bonds typical for monoterpenoids
    ]
    
    # Check if the molecule has a monoterpenoid structure and some common functional group modifications
    terpene_match = any(mol.HasSubstructMatch(pattern) for pattern in monoterpenoid_patterns)
    functional_match = any(mol.HasSubstructMatch(pattern) for pattern in functional_group_patterns)

    # Check for carbon count (should closely align with 10 carbons, slight deviations acceptable)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    if terpene_match and functional_match and (9 <= c_count <= 11):
        return True, f"Structure compatible with monoterpenoids: {c_count} carbons, matches patterns."
    
    return False, "Structure is not compatible with typical monoterpenoid patterns."