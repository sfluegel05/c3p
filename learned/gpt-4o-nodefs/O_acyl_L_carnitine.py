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
    
    # Look for the carnitine backbone pattern with specific stereochemistry
    carnitine_backbone_pattern = Chem.MolFromSmarts("[C@H](C[N+](C)(C)C)CC([O-])=O")
    if not mol.HasSubstructMatch(carnitine_backbone_pattern):
        return False, "No carnitine backbone found"
        
    # Look for the ester linkage to an acyl group allowing for variability
    # This pattern focuses on a more generic ester linkage, where R represents a variable acyl chain
    # [R] is represented by a carbon
    ester_to_acyl_pattern = Chem.MolFromSmarts("O[C@H](CC([O-])=O)OC(=O)[C]")
    if not mol.HasSubstructMatch(ester_to_acyl_pattern):
        return False, "No correct ester linkage found connecting to an acyl group"

    return True, "Contains L-carnitine backbone with ester linkage to an acyl group"