"""
Classifies: CHEBI:74716 withanolide
"""
from rdkit import Chem

def is_withanolide(smiles: str):
    """
    Determines if a molecule is a withanolide based on its SMILES string.
    Withanolides are C28 steroid lactones with modified side chains forming lactone rings and substituted derivatives.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a withanolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (None, "Invalid SMILES string")

    # Check for the steroid backbone (using a composite flexible core steroid pattern)
    steroid_patterns = [
        Chem.MolFromSmarts("C1CC(C2CCC3C4CC=CC5C4(C)CCC2C1=C35)"),  # Standard steroid scaffold
        Chem.MolFromSmarts("C1CC2CCC3C(C=CC4C3C=CC5=C4C2=C5C)C1"),  # Another steroid variation
    ]
    
    if not any(mol.HasSubstructMatch(pattern) for pattern in steroid_patterns):
        return (False, "No steroid backbone detected")
    
    # Check for lactone ring (generalized to more complex structures)
    lactone_patterns = [
        Chem.MolFromSmarts("O=C1OC=C(C)C(C)C1"),  # Variation for withanolide lactones
        Chem.MolFromSmarts("O[C@@]1([C@@H]2[C@H]3[C@@H](C(=O)O1)OC=C[C@@H]3C2)"),  # Larger lactones
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in lactone_patterns):
        return (False, "No lactone group found")
    
    # Check for common functional groups specific to withanolides
    hemiacetal_or_hydroxyl_patterns = [
        Chem.MolFromSmarts("C1(O)C(CO)CC1"),  # Pattern for hydroxyl containing rings
        Chem.MolFromSmarts("COC1=C[C@@H](O)C(O)=C(C)C1"),  # Hydroxy and ether functionalities
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in hemiacetal_or_hydroxyl_patterns):
        return (False, "Lacks characteristic oxygenated functionalities of withanolides")

    return (True, "Contains features consistent with withanolides: steroid backbone, lactone ring, and functional groups")