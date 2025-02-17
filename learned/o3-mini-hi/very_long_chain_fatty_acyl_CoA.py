"""
Classifies: CHEBI:61910 very long-chain fatty acyl-CoA
"""
"""
Classifies: Very long‐chain fatty acyl‐CoA 
Definition: A fatty acyl‐CoA in which the fatty acyl group (the acyl chain attached via a thioester bond) has a chain length greater than C22.
"""
from rdkit import Chem

def is_very_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acyl-CoA based on its SMILES string.
    The criteria are:
      1. The molecule must have a thioester group (indicating an acyl-CoA) defined by a carbonyl group immediately attached to a sulfur.
      2. The molecule must display a CoA-related moiety. Here we look for an adenine ring, a substructure in CoA.
      3. The fatty acyl chain (all carbon atoms connected to the thioester carbonyl excluding those immediately on the sulfur side)
         must have a total carbon count greater than 22.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is a very long-chain fatty acyl-CoA, False otherwise.
        str: Reason for the classification result.
    """
    # Parse the molecule from SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Step 1: Look for the thioester group.
    # We search for a carbonyl carbon attached to a sulfur.
    thioester_smarts = "[#6](=O)[S]"   # carbonyl carbon with a double bond to O and single bond to S.
    thioester = Chem.MolFromSmarts(thioester_smarts)
    if thioester is None:
        return False, "Error in thioester SMARTS"
    thioester_matches = mol.GetSubstructMatches(thioester)
    if not thioester_matches:
        return False, "No thioester group found; not an acyl-CoA"
    # Assume the first match is the thioester bond.
    thioester_match = thioester_matches[0]
    carbonyl_idx = thioester_match[0]  # the carbonyl carbon atom index
    sulfur_idx = thioester_match[1]     # the sulfur atom index

    # Step 2: Check for presence of a CoA-like moiety.
    # We search for an adenine ring, which is a core component of the CoA adenosine moiety.
    # The SMARTS pattern "c1nc2c(n1)nc[nH]2" reliably identifies an adenine substructure.
    adenine_smarts = "c1nc2c(n1)nc[nH]2"
    adenine_pat = Chem.MolFromSmarts(adenine_smarts)
    if adenine_pat is None:
        return False, "Error creating adenine SMARTS pattern"
    if not mol.HasSubstructMatch(adenine_pat):
        return False, "No CoA moiety detected (adenine fragment missing)"
        
    # Step 3: Isolate the fatty acyl chain.
    # We want to count all carbon atoms connected to the carbonyl carbon (the fatty acyl chain)
    # while ensuring that we do not traverse into the sulfur (which leads into the CoA side) or across the carbonyl oxygen.
    acyl_chain_atoms = set()
    visited = set()
    
    def dfs(atom_idx):
        visited.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        # If it is a carbon, add it to the acyl chain set.
        if atom.GetAtomicNum() == 6:
            acyl_chain_atoms.add(atom_idx)
        # Traverse neighbors.
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in visited:
                continue
            # At the starting carbonyl, restrict traversal so as not to go into the sulfur or the carbonyl oxygen.
            if atom_idx == carbonyl_idx:
                bond = mol.GetBondBetweenAtoms(atom_idx, nbr_idx)
                # Skip the sulfur attached to the carbonyl.
                if nbr_idx == sulfur_idx:
                    continue
                # If neighbor is oxygen and is double-bonded (carbonyl oxygen), skip it.
                if nbr.GetAtomicNum() == 8 and bond is not None and bond.GetBondTypeAsDouble() == 2.0:
                    continue
            # Only traverse further if the neighbor is a carbon.
            if nbr.GetAtomicNum() == 6:
                dfs(nbr_idx)
    
    # Begin DFS from the carbonyl carbon.
    dfs(carbonyl_idx)
    
    # The chain length is the number of unique carbon atoms found in the acyl fragment.
    chain_length = len(acyl_chain_atoms)
    if chain_length <= 22:
        return False, f"Fatty acyl chain has {chain_length} carbons; must be greater than 22"
    else:
        return True, f"Fatty acyl chain has {chain_length} carbons, which is >22"

# (Optional) Example usage:
# test_smiles = "CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
# result, reason = is_very_long_chain_fatty_acyl_CoA(test_smiles)
# print(result, reason)