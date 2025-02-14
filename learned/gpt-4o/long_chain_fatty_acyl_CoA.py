"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA based on its SMILES string.
    A long-chain fatty acyl-CoA results from the formal condensation of the thiol group 
    of coenzyme A with the carboxy group of any long-chain (C13 to C22) fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Enhanced Coenzyme A pattern
    coa_pattern = Chem.MolFromSmarts("C1CO[P](O)(=O)O[C@H]1C(CC(=O)NC2CCN(C(=O)S2)C(=O)NC3CCN3C(=O)[C@H](O)C(C)(C)COP(O)(=O)O)N4C=NC5=CN=CN=C54")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No coenzyme A moiety found"
    
    # Thioester linkage pattern
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"
    
    # Detect longest carbon chain for fatty acid part
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    for match in thioester_matches:
        carbon_chain = []

        # Atom index for initial thioester carbon
        atom_idx = match[0]  # carbon in thioester
        atom = mol.GetAtomWithIdx(atom_idx)
        visited = set(match)  # Start with thioester atoms
        carbon_count = 1  # Initial carbon in thioester
        max_chain_length = 0
        
        # BFS/DFS for chain traversal
        def traverse_chain(atom_idx):
            nonlocal carbon_count, max_chain_length
            queue = [atom_idx]
            while queue:
                current_idx = queue.pop()
                current_atom = mol.GetAtomWithIdx(current_idx)
                if current_atom.GetAtomicNum() == 6:  # Check for carbon
                    carbon_count += 1
                    if carbon_count > max_chain_length:
                        max_chain_length = carbon_count

                for neighbor in current_atom.GetNeighbors():
                    neighbor_idx = neighbor.GetIdx()
                    if neighbor_idx not in visited and neighbor.GetAtomicNum() == 6:
                        visited.add(neighbor_idx)
                        queue.append(neighbor_idx)

        # Start traversal from adjacent possible carbon
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                traverse_chain(neighbor.GetIdx())

        # Check if carbon chain length is within desired range
        if 13 <= max_chain_length <= 22:
            return True, f"Valid long-chain fatty acyl-CoA with {max_chain_length} carbon atoms in chain"
    
    return False, "No suitable long-chain fatty acid chain found"