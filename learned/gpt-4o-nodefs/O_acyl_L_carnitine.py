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
    
    # Look for the carnitine backbone pattern with stereo-chemistry
    carnitine_backbone_pattern = Chem.MolFromSmarts("[C@H](C[N+](C)(C)C)CC([O-])=O")
    if not mol.HasSubstructMatch(carnitine_backbone_pattern):
        return False, "No carnitine backbone found"
        
    # Look for ester linkage directly connected to an acyl group
    ester_to_acyl_pattern = Chem.MolFromSmarts("O[C@H](CC([O-])=O)C(=O)[C!H0]")  # C!H0 ensures the acyl group attaches to carbonyl carbon
    if not mol.HasSubstructMatch(ester_to_acyl_pattern):
        return False, "No correct ester linkage found connecting to an acyl group"

    return True, "Contains L-carnitine backbone with ester linkage to an acyl group"