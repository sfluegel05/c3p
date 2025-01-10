"""
Classifies: CHEBI:61905 short-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_short_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acyl-CoA based on its SMILES string.
    A short-chain fatty acyl-CoA consists of a coenzyme A moiety, a thioester linkage, and a short-chain acyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a short-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for CoA-like structure (pattern for coenzyme A moiety)
    coa_pattern = Chem.MolFromSmarts("C(=O)NCCSC")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No Coenzyme A moiety found"
    
    # Look for thioester group (C(=O)S)
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    # Count carbons attached to thioester carbon (short-chain 2 to 5 carbons)
    for match in mol.GetSubstructMatches(thioester_pattern):
        thioester_carbon_idx = match[0]
        neighbors = [n.GetIdx() for n in mol.GetAtomWithIdx(thioester_carbon_idx).GetNeighbors() if n.GetAtomicNum() == 6]
        c_count = len(neighbors)
        if 2 <= c_count <= 5:
            return True, "Contains CoA moiety and short-chain fatty acyl group"

    return False, "Does not satisfy short-chain fatty acyl criteria or invalid linkages"