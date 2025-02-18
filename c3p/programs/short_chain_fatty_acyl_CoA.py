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
    thioester_pattern = Chem.MolFromSmarts("SC(=O)C")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found" 

    # Look for Coenzyme A specific fragments with improved pattern
    # CoA includes a phosphopantetheine chain and adenine linked through a ribose phosphodiester.
    coA_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OCC1OC(n2cnc3c(ncnc23)N)C(O)C1OP(=O)([O-])O")
    if not mol.HasSubstructMatch(coA_pattern):
        return False, "Coenzyme A structure not found"

    # Check for acyl chain length from 2 to 6
    # Finds short-chain length after the thioester linkage
    for match in mol.GetSubstructMatches(thioester_pattern):
        acyl_begin = match[2]  # Carbon in SC(=O)C chain likely to be part of acyl
        chain_length = 0
        queue = [acyl_begin]
        visited = set(queue)

        while queue:
            atom_idx = queue.pop(0)
            atom = mol.GetAtomWithIdx(atom_idx)

            # Traverse chain initiated by match through carbons
            if atom.GetAtomicNum() == 6:  # Carbon
                chain_length += 1
            
            for neighbor in atom.GetNeighbors():
                n_idx = neighbor.GetIdx()
                if n_idx not in visited:
                    visited.add(n_idx)
                    queue.append(n_idx)

        # Validate short chain length
        if 2 <= chain_length <= 6:
            return True, "Contains Coenzyme A linked with short-chain fatty acid"

    return False, "No short-chain fatty acid attached to CoA"