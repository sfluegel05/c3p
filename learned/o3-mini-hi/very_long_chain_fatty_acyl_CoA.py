"""
Classifies: CHEBI:61910 very long-chain fatty acyl-CoA
"""
"""
Classifies: Very long‐chain fatty acyl‐CoA
Definition: A fatty acyl‐CoA in which the fatty acyl group (the acyl chain attached via a thioester bond)
has a chain length greater than C22.
"""
from rdkit import Chem

def is_very_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acyl-CoA based on its SMILES string.
    
    Criteria:
      1. Must have a thioester group (a carbonyl carbon attached to a sulfur).
      2. Must contain a CoA-related moiety. We check for the adenine substructure using two patterns,
         since the representation can vary.
      3. The fatty acyl chain (the chain attached to the carbonyl carbon but not including the CoA part)
         must have more than 22 carbon atoms.
    
    Args:
        smiles (str): SMILES string representing the molecule.
    
    Returns:
        bool: True if the molecule meets the criteria, False otherwise.
        str: Explanation for the classification result.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # STEP 1: Search for a thioester group.
    # Look for a pattern: a carbon (atomic number 6) double-bonded to oxygen and bonded to a sulfur.
    thioester_smarts = "[#6](=O)[S]"
    thioester = Chem.MolFromSmarts(thioester_smarts)
    if thioester is None:
        return False, "Error in thioester SMARTS definition"
    thioester_matches = mol.GetSubstructMatches(thioester)
    if not thioester_matches:
        return False, "No thioester group found; not an acyl-CoA"
    
    # Use the first found thioester substructure. In the match, the first element is the carbonyl carbon,
    # the second element is the attached sulfur.
    thioester_match = thioester_matches[0]
    carbonyl_idx = thioester_match[0]
    sulfur_idx = thioester_match[1]
    
    # STEP 2: Verify presence of a CoA moiety.
    # CoA contains an adenine group. Many SMILES of acyl-CoA include an adenine-like fragment.
    # We now try two patterns, to catch variations in atom case or substituents.
    adenine_smarts1 = "n1cnc2c(N)ncnc12"   # matches fragments like: n1cnc2c(N)ncnc12
    adenine_smarts2 = "n1cnc2ncnc12"         # alternative pattern if the extra amine is not shown
    adenine_pat1 = Chem.MolFromSmarts(adenine_smarts1)
    adenine_pat2 = Chem.MolFromSmarts(adenine_smarts2)
    if adenine_pat1 is None or adenine_pat2 is None:
        return False, "Error creating adenine SMARTS patterns"
    if not (mol.HasSubstructMatch(adenine_pat1) or mol.HasSubstructMatch(adenine_pat2)):
        return False, "No CoA moiety detected (adenine fragment missing)"
    
    # STEP 3: Isolate the fatty acyl chain.
    # We start from the carbonyl carbon and perform a DFS. We count only carbon atoms and avoid
    # traversing into the CoA part. Specifically, from the carbonyl carbon, we will not go to:
    #   - The sulfur (which leads to the CoA moiety)
    #   - The oxygen that is double-bonded (the carbonyl oxygen).
    acyl_chain_atoms = set()  # indices of carbon atoms belonging to the acyl chain
    visited = set()
    
    def dfs(atom_idx):
        visited.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        # If the atom is carbon, add it to the chain.
        if atom.GetAtomicNum() == 6:
            acyl_chain_atoms.add(atom_idx)
        # Traverse neighbors.
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in visited:
                continue
            # At the starting carbonyl carbon avoid moving toward the sulfur or the carbonyl oxygen.
            if atom_idx == carbonyl_idx:
                bond = mol.GetBondBetweenAtoms(atom_idx, nbr_idx)
                if nbr_idx == sulfur_idx:
                    continue  # avoid CoA branch via sulfur
                # Check if the neighbor is the carbonyl oxygen (double-bonded).
                if nbr.GetAtomicNum() == 8 and bond is not None and bond.GetBondTypeAsDouble() == 2.0:
                    continue
            # Only continue DFS if the neighbor is a carbon. This helps prevent entering heterocyclic rings.
            if nbr.GetAtomicNum() == 6:
                dfs(nbr_idx)
    
    dfs(carbonyl_idx)
    chain_length = len(acyl_chain_atoms)
    
    # STEP 4: Verify that the fatty acyl chain has more than 22 carbons.
    if chain_length <= 22:
        return False, f"Fatty acyl chain has {chain_length} carbons; must be greater than 22"
    else:
        return True, f"Fatty acyl chain has {chain_length} carbons, which is >22"

# (Optional) Example usage:
# test_smiles = ("CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)"
#                "[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)"
#                "[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12")
# result, reason = is_very_long_chain_fatty_acyl_CoA(test_smiles)
# print(result, reason)