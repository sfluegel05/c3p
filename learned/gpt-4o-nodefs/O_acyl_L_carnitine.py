"""
Classifies: CHEBI:75659 O-acyl-L-carnitine
"""
from rdkit import Chem

def is_O_acyl_L_carnitine(smiles: str):
    """
    Determines if a molecule is O-acyl-L-carnitine based on its SMILES string.
    O-acyl-L-carnitines have a carnitine backbone connected to an acyl group via an ester linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is O-acyl-L-carnitine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for the carnitine backbone pattern with stereochemistry
    carnitine_backbone_pattern = Chem.MolFromSmarts("[C@H](C[N+](C)(C)C)CC([O-])=O")
    if not mol.HasSubstructMatch(carnitine_backbone_pattern):
        return False, "No carnitine backbone found"
        
    # Look for ester linkage directly connected to an acyl group
    # Allow for variety in acyl chains (e.g. include flexible atom counts or types)
    ester_to_acyl_pattern = Chem.MolFromSmarts("O[C@H](CC([O-])=O)OC(=O)C")  # C after ester linkage to capture acyl part
    acyl_variability_pattern = Chem.MolFromSmarts("O=C([CX4,CX3X2])O")  # Broader pattern for ester bond
    if not mol.HasSubstructMatch(ester_to_acyl_pattern) or not mol.HasSubstructMatch(acyl_variability_pattern):
        return False, "No correct ester linkage found connecting to an acyl group"

    return True, "Contains L-carnitine backbone with ester linkage to an acyl group"