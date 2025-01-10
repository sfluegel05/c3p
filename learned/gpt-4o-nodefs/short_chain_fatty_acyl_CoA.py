"""
Classifies: CHEBI:61905 short-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_short_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acyl-CoA based on its SMILES string.
    A short-chain fatty acyl-CoA has a short fatty acid chain (<8 carbons) attached to Coenzyme A via a thioester bond.

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

    # Define updated Coenzyme A pattern
    coa_structure_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OCC1OC(C(O)C1O)N2C=NC3=C2N=CN=C3N")
    if not mol.HasSubstructMatch(coa_structure_pattern):
        return False, "No Coenzyme A structure found"

    # Look for thioester group pattern (-C(=O)S-)
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    thioester_match = mol.GetSubstructMatch(thioester_pattern)
    if not thioester_match:
        return False, "No thioester group found"

    # Find the fatty acid chain starting position just after sulfur
    fatty_acid_start_idx = thioester_match[-1]
    num_carbons = 0

    # Count carbons starting after the thioester sulfur
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            if atom.GetIdx() > fatty_acid_start_idx:
                num_carbons += 1

    if num_carbons >= 8:
        return False, f"Fatty acid chain is not short (has {num_carbons} carbons)"

    return True, "Contains short fatty acid chain connected to Coenzyme A via a thioester bond"