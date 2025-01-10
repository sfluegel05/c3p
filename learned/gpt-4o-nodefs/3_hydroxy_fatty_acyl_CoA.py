"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
from rdkit import Chem

def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES string.
    A 3-hydroxy fatty acyl-CoA has a 3-hydroxy group on the fatty acid and a CoA moiety.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if a 3-hydroxy fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for a 3-hydroxy group attached to a fatty acyl chain
    hydroxy_pattern = Chem.MolFromSmarts("[C](O)CC(=O)")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No 3-hydroxy group on fatty acyl chain found"
    
    # Simplified Coenzyme A detection focusing on key elements
    # e.g., include [SC](=O)C(=O)NCC motif as part of CoA
    coa_pattern = Chem.MolFromSmarts("C(=O)NCCSC(=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No Coenzyme A moiety found"

    return True, "Identified as 3-hydroxy fatty acyl-CoA with CoA moiety"