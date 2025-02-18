"""
Classifies: CHEBI:25409 monoterpenoid
"""
from rdkit import Chem

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.
    Monoterpenoids are derived from monoterpenes, with a C10 skeleton possibly rearranged or modified by the removal of methyl groups.

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

    # Monoterpenoid characteristics
    # Typical backbone structures including acyclic and cyclic forms
    monoterpenoid_patterns = [
        Chem.MolFromSmarts("C[C@H](C)C"),  # Basic isoprene repeat in acyclic forms
        Chem.MolFromSmarts("C1CCCCC1"),     # Cyclohexane
        Chem.MolFromSmarts("C1=CC=CC=C1"),  # Benzene rings, although rare
        Chem.MolFromSmarts("C12CCCCC1C2"),  # Bicyclic patterns
        Chem.MolFromSmarts("C1CCCC1"),      # Simple cycle variations
    ]

    # Common functional groups in monoterpenoids
    functional_group_patterns = [
        Chem.MolFromSmarts("[CX3]=[O]"),    # Carbonyl groups
        Chem.MolFromSmarts("[OX2H]"),       # Hydroxyl group
        Chem.MolFromSmarts("[NX3;H]"),      # Amine group
        Chem.MolFromSmarts("[CX3](=O)[OX2H1]"), # Acids
        Chem.MolFromSmarts("[CX3](=O)[OX2R]"),  # Esters
    ]
    
    # Check if the molecule has monoterpenoid structures and at least one functional modification
    terpene_match = any(mol.HasSubstructMatch(pattern) for pattern in monoterpenoid_patterns)
    functional_match = any(mol.HasSubstructMatch(pattern) for pattern in functional_group_patterns)

    # Checks for carbon count related to monoterpenoid core structures
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    if terpene_match and functional_match and (8 <= c_count <= 11):
        return True, f"Structure compatible with monoterpenoids: {c_count} carbons, matches patterns."
    
    return False, "Structure is not compatible with typical monoterpenoid patterns."