"""
Classifies: CHEBI:36615 triterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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
    
    # Check for triterpenoid-like patterns using SMARTS
    # Lupane, oleanane, and other base structures
    triterpenoid_patterns = [
        Chem.MolFromSmarts("C1CCCC2(C1)C3CCCCC3CC2"), # Generic pattern example
        # Placeholders for specific triterpenoid skeletons can be inserted here
        # These are not accurate representations but illustrative of concept
    ]
    
    # Check for each pattern
    for pattern in triterpenoid_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Matches triterpenoid-like pattern"
    
    # Check for common functional groups (alcohols, ketones, carboxylic acids)
    alcohol_pattern = Chem.MolFromSmarts("[CX4][OH1]")
    if mol.HasSubstructMatch(alcohol_pattern):
        return True, "Contains alcohol functional group"
    
    ketone_pattern = Chem.MolFromSmarts("C(=O)[C,C]")
    if mol.HasSubstructMatch(ketone_pattern):
        return True, "Contains ketone functional group"
    
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH1]")
    if mol.HasSubstructMatch(carboxylic_acid_pattern):
        return True, "Contains carboxylic acid functional group"
    
    return False, "Does not match typical triterpenoid structure"