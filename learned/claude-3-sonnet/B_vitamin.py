"""
Classifies: CHEBI:75769 B vitamin
"""
"""
Classifies: CHEBI:33674 B vitamin

A B vitamin is any member of the group of eight water-soluble vitamins originally thought to be a single compound 
(vitamin B) that play important roles in cell metabolism. The group comprises vitamin B1, B2, B3, B5, B6, B7, B9, 
and B12 (Around 20 other compounds were once thought to be B vitamins but are no longer classified as such).
"""

from rdkit import Chem
from rdkit.Chem import AllChem

# SMARTS patterns for common B vitamin substructures
thiazole_pattern = Chem.MolFromSmarts("c1cscn1")  # Thiazole ring (B1)
isoalloxazine_pattern = Chem.MolFromSmarts("c1nc2c(nc1N)c(=O)nc(=O)n2")  # Isoalloxazine ring (B2)
pyridine_pattern = Chem.MolFromSmarts("c1ncccc1")  # Pyridine ring (B3, B5, B6)
pterin_pattern = Chem.MolFromSmarts("c1nc2c(nc(N)nc2n1)N")  # Pterin ring (B9)
corrin_pattern = Chem.MolFromSmarts("[Co]1234(N5C6=C(C)C7=[N+]1C(=CC1=[N+]2C(=C(C)C2=[N+]3[C@]5([H])[C@@](C)(CC(N)=O)[C@@H]4CCC(N)=O)[C@@](C)(CC(N)=O)[C@@H]3CCC(N)=O)C(C)(C)[C@@H]1CCC(N)=O")  # Corrin ring (B12)

def is_B_vitamin(smiles: str):
    """
    Determines if a molecule is a B vitamin based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a B vitamin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for common B vitamin substructures
    if mol.HasSubstructMatch(thiazole_pattern):
        return True, "Contains thiazole ring (B1)"
    if mol.HasSubstructMatch(isoalloxazine_pattern):
        return True, "Contains isoalloxazine ring (B2)"
    if mol.HasSubstructMatch(pyridine_pattern):
        # Further checks for B3, B5, B6
        pass
    if mol.HasSubstructMatch(pterin_pattern):
        return True, "Contains pterin ring (B9)"
    if mol.HasSubstructMatch(corrin_pattern):
        return True, "Contains corrin ring (B12)"
    
    # If no common substructures found, check for other B vitamin patterns
    # ...
    
    return False, "No B vitamin substructures found"