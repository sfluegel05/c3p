"""
Classifies: CHEBI:61905 short-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_short_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acyl-CoA based on its SMILES string.

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
    
    # Look for thioester linkage: S-C(=O)
    thioester_pattern = Chem.MolFromSmarts("S-C(=O)")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"
    
    # Look for Coenzyme A specific fragments (pattern recognition)
    coA_pattern = Chem.MolFromSmarts("O=C(O)CCNC(=O)C[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)O")  # Simplified CoA part
    if not mol.HasSubstructMatch(coA_pattern):
        return False, "Coenzyme A structure not found"
    
    # Check for short-chain acyl group (e.g., 2 to 6 carbons in acyl chain)
    acyl_chain_pattern = Chem.MolFromSmarts("CC(=O)S")
    matches = mol.GetSubstructMatches(acyl_chain_pattern)
    
    for match in matches:
        # Check the length of carbon chain directly attached to the S=C(=O) group
        acyl_chain_len = 0
        atom_indices = [atom_idx for atom_idx in match]
        sulfur_idx, carbon_idx, o_carbon_idx = atom_indices[0], atom_indices[1], atom_indices[2]
        
        for bond in mol.GetBondWithIdx(mol.GetBondBetweenAtoms(sulfur_idx, carbon_idx).GetIdx()).GetConnectedAtoms():
            if bond.GetAtomicNum() == 6:
                acyl_chain_len += 1
        
        if 2 <= acyl_chain_len <= 6:
            return True, "Contains Coenzyme A with short-chain fatty acid"
    
    return False, "No short-chain fatty acid attached to CoA"