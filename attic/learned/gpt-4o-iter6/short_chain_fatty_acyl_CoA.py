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
        bool: True if the molecule is a short-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify the Coenzyme A pattern
    # Adjusting the pattern to allow flexibility in the `C` and `O` groups that accompany CoA
    coa_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)[C@H](O)C")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No identifiable Coenzyme A moiety found"

    # Detect thioester and short chain (2-5 carbons)
    thioester_pattern = Chem.MolFromSmarts("C(=O)SC")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)

    for match in thioester_matches:
        thioester_carbon = match[0]  # Thioester carbon
        visited = set()
        carbon_chain = []
        atom_queue = [(thioester_carbon)]
        
        while atom_queue:
            atom_idx = atom_queue.pop()
            if atom_idx in visited:
                continue
            visited.add(atom_idx)

            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 6:  # Carbon
                carbon_chain.append(atom_idx)
                if len(carbon_chain) > 5:  # Exceeded short-chain length
                    break
                for neighbor in atom.GetNeighbors():
                    neighbor_idx = neighbor.GetIdx()
                    if neighbor_idx not in visited:  # Check all neighbors
                        atom_queue.append(neighbor_idx)

        if 2 <= len(carbon_chain) <= 5:  # Short-chain (2 to 5 carbons)
            return True, "CoA moiety and short-chain fatty acyl thioester linkage present"

    return False, "Does not satisfy short-chain fatty acyl criteria, or invalid linkages"