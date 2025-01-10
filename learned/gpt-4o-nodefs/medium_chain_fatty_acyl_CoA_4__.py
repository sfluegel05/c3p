"""
Classifies: CHEBI:90546 medium-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem

def is_medium_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA(4-) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a medium-chain fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Improved Coenzyme A SMARTS pattern (more comprehensive)
    # Includes ribose, phosphate groups, and important functionalities
    coa_comprehensive_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)CCSC(=O)CNC(=O)[C@H]1[C@@H](O)[C@@H](O)[C@H](COP(O)(=O)OP(O)(=O)OC[C@H]2O[C@@H](n3cnc4c(N)ncnc34)[C@@H]2O)O1")
    if not mol.HasSubstructMatch(coa_comprehensive_pattern):
        return False, "Coenzyme A moiety not found"

    # Check for thioester linkage (C(=O)S)
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    # Check for flexible medium acyl chain
    # Using custom logic to count carbon atoms in a chain-like structure
    fatty_acid_pattern = Chem.MolFromSmarts("C(=O)SCC")
    matches = mol.GetSubstructMatches(fatty_acid_pattern)
    for match in matches:
        atoms = [mol.GetAtomWithIdx(idx) for idx in match]
        chain_length = 0
        atom_queue = [atoms[-1]]  # Start from the S-C-C pattern
        visited = set(match)
        while atom_queue:
            atom = atom_queue.pop()
            chain_length += 1 if atom.GetAtomicNum() == 6 else 0
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in visited and neighbor.GetAtomicNum() in (1, 6):  # Allow H and C only
                    visited.add(neighbor.GetIdx())
                    atom_queue.append(neighbor)
        if 6 <= chain_length <= 12:
            return True, "Contains medium-chain fatty acyl linked to Coenzyme A"
    
    return False, "No medium-chain fatty acyl chain found"