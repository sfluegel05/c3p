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
    The criteria are:
      1. The molecule must have a thioester group (a carbonyl carbon attached to a sulfur atom).
      2. The molecule must display a CoA-related moiety. We look for an adenine substructure,
         but since the representation varies, we try two slightly different SMARTS patterns.
      3. The fatty acyl chain (the chain attached to the thioester's carbonyl carbon, 
         but not including the sulfur side leading into CoA) must have a carbon count greater than 22.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule meets the very long-chain fatty acyl-CoA criteria, False otherwise.
        str: Reason for the classification result.
    """
    # Parse molecule from SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1: Look for a thioester group.
    # We search for a carbonyl carbon immediately attached to a sulfur.
    thioester_smarts = "[#6](=O)[S]"  # carbonyl group attached to sulfur
    thioester = Chem.MolFromSmarts(thioester_smarts)
    if thioester is None:
        return False, "Error in thioester SMARTS definition"
    thioester_matches = mol.GetSubstructMatches(thioester)
    if not thioester_matches:
        return False, "No thioester group found; not an acyl-CoA"
    
    # Use the first thioester match.
    thioester_match = thioester_matches[0]
    carbonyl_idx = thioester_match[0]  # this is the carbonyl carbon index
    sulfur_idx = thioester_match[1]     # this is the sulfur index
    
    # Step 2: Check for the presence of a CoA-related moiety.
    # Coenzyme A contains an adenine substructure. However, the representation in SMILES can vary,
    # so we test two patterns and accept if either is present.
    adenine_smarts1 = "n1cnc2ncnc12"         # pattern using all lowercase atoms
    adenine_smarts2 = "c1nc2c(n1)nc[nH]2"      # alternative pattern
    adenine_pat1 = Chem.MolFromSmarts(adenine_smarts1)
    adenine_pat2 = Chem.MolFromSmarts(adenine_smarts2)
    if adenine_pat1 is None or adenine_pat2 is None:
        return False, "Error creating adenine SMARTS patterns"
    if not (mol.HasSubstructMatch(adenine_pat1) or mol.HasSubstructMatch(adenine_pat2)):
        return False, "No CoA moiety detected (adenine fragment missing)"
    
    # Step 3: Isolate and count the fatty acyl chain.
    # We perform a DFS starting from the carbonyl carbon, counting only carbon atoms 
    # but ensuring we do not traverse into the sulfur (which would lead into the CoA part)
    # or across the carbonyl oxygen.
    acyl_chain_atoms = set()
    visited = set()
    
    def dfs(atom_idx):
        visited.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        # If the atom is a carbon, record it as part of the acyl chain.
        if atom.GetAtomicNum() == 6:
            acyl_chain_atoms.add(atom_idx)
        # Traverse neighboring atoms.
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in visited:
                continue
            # At the starting carbonyl carbon, we do not go toward the sulfur or the carbonyl oxygen.
            if atom_idx == carbonyl_idx:
                bond = mol.GetBondBetweenAtoms(atom_idx, nbr_idx)
                # Skip the sulfur atom directly attached.
                if nbr_idx == sulfur_idx:
                    continue
                # Skip the double-bonded oxygen (the carbonyl oxygen).
                if nbr.GetAtomicNum() == 8 and bond is not None and bond.GetBondTypeAsDouble() == 2.0:
                    continue
            # Continue the DFS only if the neighbor is carbon.
            if nbr.GetAtomicNum() == 6:
                dfs(nbr_idx)
                
    dfs(carbonyl_idx)
    
    chain_length = len(acyl_chain_atoms)
    if chain_length <= 22:
        return False, f"Fatty acyl chain has {chain_length} carbons; must be greater than 22"
    else:
        return True, f"Fatty acyl chain has {chain_length} carbons, which is >22"

# (Optional) Example usage:
# test_smiles = "CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
# result, reason = is_very_long_chain_fatty_acyl_CoA(test_smiles)
# print(result, reason)