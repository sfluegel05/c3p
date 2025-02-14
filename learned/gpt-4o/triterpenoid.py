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
    
    # Define skeletal structures typical of triterpenoids
    triterpenoid_patterns = [
        Chem.MolFromSmarts("C1CCC2CC[C@H]3[C@]4(C)CC[C@@H]5C4CC[C@@H]3[C@@H]2C1"),
        Chem.MolFromSmarts("C1CC2(CC[C@@H]3CC4(C)CCC5C4CCC3C2CCC1C5)C"),
        Chem.MolFromSmarts("C1=CC2=C3C=CC4=C2C1C5=C3C(C=C6[C@H](C5)CC(=CC6)C)C4"),
    ]
    
    # Check if matches any typical triterpenoid patterns
    for pattern in triterpenoid_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Matches typical triterpenoid skeletal pattern"
    
    # Check for common triterpenoid functional group patterns within specific contexts
    specific_alcohol_pattern = Chem.MolFromSmarts("[C&!R][C@H](O)C")
    if mol.HasSubstructMatch(specific_alcohol_pattern):
        return True, "Contains alcohol functional group in sterane ring"

    specific_ketone_pattern = Chem.MolFromSmarts("[C&R2](=O)[C&R2]")
    if mol.HasSubstructMatch(specific_ketone_pattern):
        return True, "Contains ketone functional group in sterane ring"

    specific_carboxylic_acid_pattern = Chem.MolFromSmarts("[C&R2](=O)[OH]")
    if mol.HasSubstructMatch(specific_carboxylic_acid_pattern):
        return True, "Contains carboxylic acid functional group in sterane ring"

    return False, "Does not match typical triterpenoid structure"