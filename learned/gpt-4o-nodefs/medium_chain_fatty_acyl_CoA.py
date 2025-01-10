"""
Classifies: CHEBI:61907 medium-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_medium_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a medium-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Coenzyme A moiety pattern
    coa_pattern = Chem.MolFromSmarts("OP(O)(=O)OCC1OC(COP(O)(=O)O)C1OP(O)(=O)O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety found"

    # Thioester bond pattern
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    # Finding Acyl Chain Length
    connected_parts = thioester_pattern.GetAtoms()
    for start_atom in connected_parts:
        # Consider the connectivity from thioester justification
        chain_length = 0
        for neighbor in start_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon atoms
                chain_length += 1

    # Check if there is any medium chain (6-12 carbons) acyl component
    if 6 <= chain_length <= 12:
        return True, "Contains a CoA moiety with a medium-length acyl chain"

    return False, "Fatty acyl chain not of medium length"