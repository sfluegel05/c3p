"""
Classifies: CHEBI:61655 steroid saponin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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

    # Modified steroid backbone pattern; more versatile to account for known variations
    steroid_patterns = [
        Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CC4=CC(=O)CC=C4C3"),
        Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2C=C4C=CC(=O)CC4C3"),  # Provide different steroid skeletons
        # Add more versatile steroid patterns if necessary
    ]
    if not any(mol.HasSubstructMatch(p) for p in steroid_patterns):
        return False, "No steroid backbone found"
        
    # Enhanced glycosidic linkage pattern allowing for variable sugar attachments
    # This pattern takes into account the complexity and variability of glycosidic linkages
    glycosidic_linkage = Chem.MolFromSmarts("[C,O]-[O]-[C]")
    if not mol.HasSubstructMatch(glycosidic_linkage):
        return False, "No glycosidic linkage found"
        
    # Improved pattern for sugar moieties; aims to capture general sugar-like structures
    sugar_patterns = [
        Chem.MolFromSmarts("C(O)C(O)C(O)C"),  # Simple pattern for sugars
        Chem.MolFromSmarts("C(C(CO)O)O"),      # Patterns for common sugar variants
        # Additional sugar patterns could be defined here
    ]
    if not any(mol.HasSubstructMatch(p) for p in sugar_patterns):
        return False, "No sugar moiety found"

    # Additional checks for molecular weight and ring system
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight < 400:  # Steroid saponins typically have higher molecular weights
        return False, "Molecular weight too low for steroid saponin"

    return True, "Contains steroid backbone with glycosidic linkage to sugar moieties"