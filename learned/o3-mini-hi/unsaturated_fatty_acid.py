"""
Classifies: CHEBI:27208 unsaturated fatty acid
"""
"""
Classifies: Unsaturated Fatty Acid
Definition:
    An unsaturated fatty acid is defined herein as a molecule that meets all of the following:
      1. The molecule is valid.
      2. It contains exactly one carboxylic acid group (COOH or COO-), whose carbonyl carbon is attached 
         to exactly one carbon (i.e. terminal acid functionality).
      3. Starting from that sole chain‐attached carbon, the longest simple path (via carbon–carbon bonds) 
         is at least 6 carbons in length (this heuristically ensures a long aliphatic chain).
      4. The molecule contains at least one carbon–carbon double bond (C=C) and no triple bonds.
      
Note:
  – This is a heuristic approach. Some free fatty acids may contain rings or substituents that complicate 
    chain detection. Here we allow rings in parts of the molecule (so that e.g. leukotriene A4, juvenile hormone I acid,
    etc. are not rejected solely for having rings) but we require that the acid group is terminal and that one can trace 
    a sufficiently long carbon‐only path from it.
      
Examples that should be classified as unsaturated fatty acids include many with long aliphatic chains plus a C=C bond, 
while molecules having extra acid groups, triple bonds, or an ambiguous/branched acid attachment are rejected.
"""

from rdkit import Chem

def is_unsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acid based on its SMILES string.
    
    The criteria are:
      1. Molecule must be valid.
      2. It must contain exactly one carboxylic acid group (protonated or deprotonated) as defined by SMARTS.
      3. The carboxyl carbon (the one in C(=O)[O;H1,O-]) should be terminal,
         meaning it is attached to exactly one carbon atom.
      4. The alkyl chain (as measured by the longest simple path of carbon atoms starting from
         the acid’s neighbor) is at least 6 carbons long.
      5. The molecule must contain at least one C=C double bond and must not contain any C#C triple bonds.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as an unsaturated fatty acid, False otherwise.
        str: Explanation of the result.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a carboxylic acid group.
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1,O-]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "Does not contain a carboxylic acid group"
    if len(acid_matches) != 1:
        return False, "Molecule should have exactly one carboxylic acid group"
    
    # In our SMARTS, atom 0 is the carbonyl carbon.
    acid_match = acid_matches[0]
    acid_c_idx = acid_match[0]
    acid_c = mol.GetAtomWithIdx(acid_c_idx)
    
    # For a terminal acid, the acid carbon should be attached to exactly one carbon.
    carbon_neighbors = [nbr.GetIdx() for nbr in acid_c.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1:
        return False, "The carboxyl group is not terminal (acid carbon is attached to more than one carbon)"
    chain_start_idx = carbon_neighbors[0]
    
    # Now, we try to measure the “fatty acid chain length”:
    # We define it as the length (number of carbons) of the longest simple path (no revisiting atoms)
    # starting from the acid group’s attached carbon. We only consider atoms with atomic number 6.
    # (This DFS-based approach is heuristic.)
    def longest_carbon_path(current_idx, visited):
        visited.add(current_idx)
        current_atom = mol.GetAtomWithIdx(current_idx)
        max_length = 1  # count the current carbon
        for nbr in current_atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            nbr_idx = nbr.GetIdx()
            if nbr_idx in visited:
                continue
            # Recurse to get path length from neighbor.
            new_length = 1 + longest_carbon_path(nbr_idx, visited.copy())
            if new_length > max_length:
                max_length = new_length
        return max_length

    chain_length = longest_carbon_path(chain_start_idx, set())
    if chain_length < 6:
        return False, f"Alkyl chain too short for a fatty acid (chain length = {chain_length})"
    
    # Check for unsaturation:
    # Require at least one C=C double bond.
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(double_bond_pattern):
        return False, "Contains no carbon–carbon double bonds (no unsaturation)"
    
    # Reject any triple bonds.
    triple_bond_pattern = Chem.MolFromSmarts("C#C")
    if mol.HasSubstructMatch(triple_bond_pattern):
        return False, "Contains carbon–carbon triple bonds which are not allowed"
    
    return True, ("Contains a terminal carboxyl group, a sufficiently long carbon chain (length = "
                  f"{chain_length}), and at least one C=C double bond (unsaturation)")

# Example usage (you can test with some of the examples listed):
if __name__ == '__main__':
    test_smiles = [
        # True positives
        "C(C(O)=O)C[C@@H](/C=C/C=C\\C/C=C\\C/C=C\\C=C\\[C@H](C/C=C\\CC)O)O",  # resolvin D6
        "CCCCCCC/C=C/CCCCCCC/C(=O)O",  # simplified (2E,4E,6E,8E,10E)-octadecapentaenoic acid-like
        "CC(C)=CCCC(C)=CCCC(C)=CC(O)=O",  # farnesoic acid-like
        "OC(=O)CCC=C",  # pent-4-enoic acid (should be rejected due to chain length)
        "[H]C(C)=CC(O)=O",  # isocrotonic acid-like (chain too short)
        "OC(=O)CC#CCCCCC#CCCCCCC",  # molecule with triple bonds (should be rejected)
    ]
    
    for smi in test_smiles:
        result, reason = is_unsaturated_fatty_acid(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n{'-'*50}")