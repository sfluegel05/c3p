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

    # Terpene-derived core patterns (cyclic variations)
    terpene_patterns = [
        Chem.MolFromSmarts("C1=CC=CC(C1)C(C)C"),  # Linear monoterpenes with rearrangement
        Chem.MolFromSmarts("C1=CCCC=C1C(C)C"),    # Cyclohexene-based variations
        Chem.MolFromSmarts("C12CCCCC1C2C(C)C"),   # Bicyclic monoterpenoids
    ]
    
    # Functional groups typically found in monoterpenoids
    functional_group_patterns = [
        Chem.MolFromSmarts("[CX3]=[O]"), # Carbonyl
        Chem.MolFromSmarts("[CX4]O"),    # Alcohol
        Chem.MolFromSmarts("[CX3](=O)[OX2]"), # Ester-like moiety
        Chem.MolFromSmarts("[OX2][CX3]=[OX1]"), # Carboxylic acid
    ]

    # Verify carbon count, focusing on the heuristic range for monoterpenoids
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Primary logic checking for motifs and ensuring flexibility
    terpene_match = any(mol.HasSubstructMatch(tp) for tp in terpene_patterns)
    functional_match = any(mol.HasSubstructMatch(fg) for fg in functional_group_patterns)
    
    if terpene_match and functional_match and 9 <= c_count <= 11:
        return True, f"Structure compatible with monoterpenoids: {c_count} carbons, observed structure"
    
    return False, "Structure is not compatible with typical monoterpenoid patterns."