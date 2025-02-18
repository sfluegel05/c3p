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

    # Look for Coenzyme A specific fragments
    # CoA includes a phosphopantetheine chain and adenine linked through a ribose phosphodiester
    coA_pattern = Chem.MolFromSmarts("OP(=O)(O)OC[C@H]1O[C@H]([C@@H](O)[C@H]1OP(=O)(O)O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coA_pattern):
        return False, "Coenzyme A structure not found"

    # Look for thioester linkage: S-C(=O) for fatty acyl attachment
    thioester_pattern = Chem.MolFromSmarts("SC(=O)C")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester linkage found"

    # Check for acyl chain length from 2 to 6 carbons
    for match in thioester_matches:
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
            return True, "Contains Coenzyme A linked with a short-chain fatty acid"

    return False, "No short-chain fatty acid attached to CoA"