"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
from rdkit import Chem

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
    
    # Coenzyme A-specific SMARTS pattern: Recognizes phosphorylated ribose, adenosine, phosphopantetheine
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(O)(=O)OC[C@H]1O[C@H](CO[P](O)(O)=O)[C@H](O)[C@H]1OP(O)(O)=O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No coenzyme A moiety found"
    
    # Thioester linkage pattern for recognizing fatty acyl attachment
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCCN")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester linkage found"

    # Detect longest carbon chain attached to thioester
    for match in thioester_matches:
        carbon_chain = []

        # Atom index for thioester carbon
        atom_idx = match[0]
        visited = set(match)
        carbon_count = 0
        
        def traverse_chain(atom_idx):
            nonlocal carbon_count
            queue = [atom_idx]
            while queue:
                current_idx = queue.pop()
                current_atom = mol.GetAtomWithIdx(current_idx)
                visited.add(current_idx)
                
                if current_atom.GetAtomicNum() == 6:  # Check for carbon
                    carbon_chain.append(current_idx)
                    carbon_count += 1

                for neighbor in current_atom.GetNeighbors():
                    neighbor_idx = neighbor.GetIdx()
                    if neighbor_idx not in visited and neighbor.GetAtomicNum() == 6:
                        queue.append(neighbor_idx)

        for neighbor in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                traverse_chain(neighbor.GetIdx())

        # Check carbon chain length
        if 13 <= carbon_count <= 22:
            return True, f"Valid long-chain fatty acyl-CoA with {carbon_count} carbon atoms in chain"
    
    return False, "No suitable long-chain fatty acid chain found"