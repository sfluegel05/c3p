"""
Classifies: CHEBI:36615 triterpenoid
"""
from rdkit import Chem

def is_triterpenoid(smiles: str):
    """
    Determines if a molecule is a triterpenoid based on its SMILES string.
    Triterpenoids are terpenoids derived from a triterpene with a varied 
    C30 carbon skeleton that might be rearranged or modified.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for at least 30 carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 30:
        return False, f"Contains only {c_count} carbon atoms, requires at least 30 for triterpenoid"
    
    # Define skeletal structures typical of triterpenoids, focusing on sterane-like structures
    triterpenoid_patterns = [
        Chem.MolFromSmarts("[C@@]1(C[C@@H]2CC[C@]3(C)[C@]4(CC[C@@H](O)[C@]5(C)C)C=C[C@@H]3[C@H]4CC[C@@]25C)C"),
        Chem.MolFromSmarts("C1CCC2(C(C)CCC3[C@H]4CC[C@@H]5[C@@](C4CCC3)[C@H](O)[C@]5(C)CO)C2C1"), 
        # Additionally decorated ring structures
    ]
    
    # Check against complex triterpenoid skeletal patterns
    for pattern in triterpenoid_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Matches triterpenoid skeletal pattern"
    
    # Check for common triterpenoid functional group patterns within specific contexts
    common_functional_groups = [
        Chem.MolFromSmarts("[CX3](=O)[OX2H1]"), # Carboxylic acid
        Chem.MolFromSmarts("[CX3](=O)[C]"),    # Ketone
        Chem.MolFromSmarts("[CX4][OX2H]")      # Alcohol
    ]
    
    matches = sum(mol.HasSubstructMatch(fg) for fg in common_functional_groups)
    if matches > 0:
        return True, f"Contains {matches} triterpenoid-like functional groups"
    
    return False, "Does not match typical triterpenoid structure"