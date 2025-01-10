"""
Classifies: CHEBI:61910 very long-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_very_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acyl-CoA based on its SMILES string.
    A very long-chain fatty acyl-CoA typically features an acyl chain with at least 22 carbons
    attached to a CoA moiety via a thioester bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a very long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify the thioester group pattern
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester group found"

    # Full CoA moiety SMARTS pattern
    # Including the adenylate as part of CoA detection is vital
    # This includes key features of ribose, phosphate linkage, and adenine
    coa_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)P(O)(O)=O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not detected"

    # Count carbon atoms in the longest aliphatic chain part of thioester
    longest_chain_length = 0
    for frag in mol.GetSubstructMatches(thioester_pattern):
        starting_atom = frag[0]  # carbonyl carbon
        atoms_visited = {starting_atom}
        carbons_in_chain = set()
        atom_stack = [mol.GetAtomWithIdx(starting_atom)]
        
        while atom_stack:
            atom = atom_stack.pop()
            if atom.GetAtomicNum() == 6:
                carbons_in_chain.add(atom.GetIdx())

            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in atoms_visited:
                    atoms_visited.add(neighbor.GetIdx())
                    atom_stack.append(neighbor)
        
        chain_length = len(carbons_in_chain)
        longest_chain_length = max(longest_chain_length, chain_length)
    
    if longest_chain_length < 22:
        return False, "Aliphatic chain not sufficiently long (very long-chain fatty acids typically have ≥22 carbons)"

    return True, "Contains a long fatty acyl chain and CoA moiety with thioester bond"

# Example test
example_smiles = "CCCCCCCCCCCCCCCCCCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)P(O)(O)=O)n1cnc2c(N)ncnc12"
result, reason = is_very_long_chain_fatty_acyl_CoA(example_smiles)
print(result, reason)