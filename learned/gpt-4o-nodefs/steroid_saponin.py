"""
Classifies: CHEBI:61655 steroid saponin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_steroid_saponin(smiles: str):
    """
    Determines if a molecule is a steroid saponin based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a steroid saponin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Better patterns for the steroid backbone
    steroid_patterns = [
        Chem.MolFromSmarts("C1CCC2C(C1)CCC3C(C2)CCC4)(C=C3)CCC4"),  # Cholestane skeleton
        Chem.MolFromSmarts("C1CCC2C(C1)C3CCC4(C2)C=CC(=O)CC4O"),    # Spirostane skeleton
        Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4CCC3OC4"),
        # More spirostane variations can be added
    ]
    if not any(mol.HasSubstructMatch(p) for p in steroid_patterns):
        return False, "No suitable steroid backbone found"
        
    # Enhanced pattern for glycosidic linkage
    glycosidic_linkage = Chem.MolFromSmarts("[O]-[C]([O,CO])[O]")
    if not mol.HasSubstructMatch(glycosidic_linkage):
        return False, "No glycosidic linkage found"
        
    # Improved patterns for sugar moieties with variability in the structures
    sugar_patterns = [
        Chem.MolFromSmarts("C(O[C,N])[C@H](O)[C@H](O)[C@H](O)"),  # Common sugar linkages
        Chem.MolFromSmarts("O[C@H]1[C@@H](O)[C@H](O[C@@H]2O)[C@@H](O)[C@H](O)[C@@H]1O"),  # Cyclic sugars
        # Consider adding patterns for various substituents like acetyl groups if needed
    ]
    if not any(mol.HasSubstructMatch(p) for p in sugar_patterns):
        return False, "No sugar moiety found"

    # Additional checks: the presence of multiple sugar moieties connected to the steroid core
    # counting the number of times glycosidic linkage appears might indicate multiple sugars
    if mol.GetSubstructMatches(glycosidic_linkage) < 2:
        return False, "Insufficient glycosidic linkages for steroid saponin (usually has more than one sugar)"

    return True, "Contains a suitable steroid backbone linked with sugar moieties via glycosidic bonds"