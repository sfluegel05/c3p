"""
Classifies: CHEBI:75659 O-acyl-L-carnitine
"""
from rdkit import Chem

def is_O_acyl_L_carnitine(smiles: str):
    """
    Determines if a molecule is O-acyl-L-carnitine based on its SMILES string.

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

    # Look for carnitine backbone pattern
    carnitine_pattern = Chem.MolFromSmarts("C[N+](C)(C)C")
    if not mol.HasSubstructMatch(carnitine_pattern):
        return False, "No carnitine backbone found"
        
    # Look for ester linkage connecting carnitine
    ester_linkage_pattern = Chem.MolFromSmarts("O[C@H](CC([O-])=O)C[N+](C)(C)C")
    if not mol.HasSubstructMatch(ester_linkage_pattern):
        return False, "No ester linkage found between carnitine and acyl group"

    # Must have at least one acyl group
    acyl_group_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(acyl_group_pattern):
        return False, "No acyl group present"

    return True, "Contains carnitine backbone with ester linkage to an acyl group"