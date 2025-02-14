"""
Classifies: CHEBI:61907 medium-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_medium_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA based on its SMILES string.
    A medium-chain fatty acyl-CoA is a fatty acyl-CoA with a fatty acid chain of 6-12 C atoms.

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
    
    # Comprehensive pattern for Coenzyme A core structure
    coa_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)[C@H](O)[C@H](COP(=O)(O)OP(=O)(O)OC[C@@H]1O[C@@H](O)[C@@H](O)C1O)n1cnc2c(ncnc2n1)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing Coenzyme A structure"

    # Thioester linkage pattern part of fatty acid-CoA connection
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester linkage found"

    # Analyze carbon chain length attached to thioester linkage
    # Assume generally linear chain for analysis; adjustment needed for complex branching
    for match in thioester_matches:
        thioester_atom_index = match[1]  # S in C(=O)S
        carbon_chain = []
        atom_queue = [atom.GetIdx() for atom in mol.GetAtomWithIdx(thioester_atom_index).GetNeighbors() if atom.GetAtomicNum() == 6]
        visited = set(atom_queue)
        while atom_queue:
            atom_idx = atom_queue.pop(0)
            carbon_chain.append(atom_idx)
            for neighbor in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                    visited.add(neighbor.GetIdx())
                    atom_queue.append(neighbor.GetIdx())
        
        chain_length = len(carbon_chain)
        if 6 <= chain_length <= 12:
            return True, "Contains medium-chain fatty acyl-CoA structure with correct components and chain length"
    
    return False, "Fatty acid chain length not within 6-12 carbons or incorrect"