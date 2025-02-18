"""
Classifies: CHEBI:27208 unsaturated fatty acid
"""
"""
Classifies: Unsaturated Fatty Acid
Definition:
    A fatty acid is assumed to have a terminal carboxylic acid group attached to a long 
    aliphatic chain. An unsaturated fatty acid must also contain at least one carbon–carbon 
    double (C=C) or triple (C#C) bond. Additionally, we require that the molecule is “simple” 
    (acyclic) so that it represents a free fatty acid rather than a complex conjugated structure.
    
Improvements over the previous attempt:
  - Reject molecules with rings.
  - Require exactly one (terminal) carboxylic acid group.
  - Ensure that the acid group is terminal by checking that its carbonyl carbon is connected 
    to exactly one aliphatic carbon.
  - Traverse that carbon chain and require a minimum chain length.
  - Still require at least one unsaturation.
  
Note: This is a heuristic approach. Some borderline cases may still be mis‐classified.
"""

from rdkit import Chem

def is_unsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acid based on its SMILES string.
    
    The criteria are:
      1. The molecule must be valid and acyclic (no rings).
      2. It must contain exactly one carboxylic acid group (COOH or COO- form).
      3. The acid group must be terminal (attached to a long, predominantly carbon chain).
      4. The chain (attached to the acid) should be at least 6 carbons long.
      5. The molecule must contain at least one C=C or C#C bond outside of the COOH.
      
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as an unsaturated fatty acid; False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Reject molecules with rings (most free fatty acids are acyclic)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings; expected a linear fatty acid chain"
    
    # Define a SMARTS pattern for a carboxylic acid group.
    # This pattern will match both protonated (COOH) and deprotonated (COO-) forms.
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1,O-]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "Does not contain a carboxylic acid group"
    
    # For simplicity, require exactly one acid group
    if len(acid_matches) != 1:
        return False, "Molecule should have exactly one carboxylic acid group"
    
    # In our SMARTS the first atom is the carbonyl carbon.
    acid_match = acid_matches[0]
    acid_c_idx = acid_match[0]
    acid_c = mol.GetAtomWithIdx(acid_c_idx)
    
    # Find the non-oxygen neighbor(s) of this acid carbon.
    # The acid carbon should be attached to one oxygen (in the acid group) and one carbon (the chain)
    carbon_neighbors = [nbr.GetIdx() for nbr in acid_c.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1:
        return False, "The carboxyl group is not terminal (expected one alkyl chain attachment)"
    chain_start_idx = carbon_neighbors[0]
    
    # Traverse the chain (only through carbon atoms) using DFS to get the longest connected chain length.
    def dfs(current_idx, visited):
        visited.add(current_idx)
        current_atom = mol.GetAtomWithIdx(current_idx)
        max_length = 1  # count current carbon
        for nbr in current_atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:  # only continue through carbon atoms
                continue
            nbr_idx = nbr.GetIdx()
            if nbr_idx in visited:
                continue
            # Continue traversal
            chain_length = 1 + dfs(nbr_idx, visited.copy())
            if chain_length > max_length:
                max_length = chain_length
        return max_length

    chain_length = dfs(chain_start_idx, set())
    if chain_length < 6:
        return False, f"Alkyl chain too short for a fatty acid (chain length = {chain_length})"
    
    # Check for unsaturation (C=C or C#C). We expect at least one in the molecule.
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    triple_bond_pattern = Chem.MolFromSmarts("C#C")
    has_double = mol.HasSubstructMatch(double_bond_pattern)
    has_triple = mol.HasSubstructMatch(triple_bond_pattern)
    
    if not (has_double or has_triple):
        return False, "Contains no carbon–carbon double or triple bonds (no unsaturation)"
    
    return True, "Contains a terminal carboxyl group, a long aliphatic chain, and at least one C=C or C#C bond (unsaturation)"

# Example usage:
if __name__ == '__main__':
    # Try a few examples; these are some known fatty acids.
    test_smiles = [
        "OC(=O)CCC=C",  # pent-4-enoic acid (valid unsaturated fatty acid)
        "CCCCCCCCCCCCCCCC(=O)O",  # long saturated fatty acid
        "C1CCCC1C(=O)O",  # cyclic acid (should be rejected)
        "CCCC(=O)O",  # too short chain
        "CC(C)=CCCC(C)=CCCC(C)=CC(O)=O",  # farnesoic acid (has unsaturation)
    ]
    
    for smi in test_smiles:
        result, reason = is_unsaturated_fatty_acid(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n{'-'*50}")