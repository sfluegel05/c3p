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
      1. The molecule must have a thioester group (indicating a fatty acyl-CoA).
      2. The molecule must display a CoA-related moiety (here we search for part of the adenosine ring).
      3. The fatty acyl chain (all carbon atoms connected to the carbonyl carbon, excluding the sulfur and carbonyl oxygen) must have a total of more than 22 carbon atoms.
      
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a very long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse the molecule from SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Step 1: Look for the thioester group.
    # We search for a pattern "[#6](=O)[S]" meaning a carbonyl carbon attached to a sulfur.
    thioester_smarts = "[#6](=O)[S]"
    thioester = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester)
    if not thioester_matches:
        return False, "No thioester group found; not an acyl-CoA"
    # Assume the first match is the thioester bond.
    thioester_match = thioester_matches[0]
    carbonyl_idx = thioester_match[0]  # the carbonyl carbon atom index
    sulfur_idx = thioester_match[1]     # the sulfur atom index

    # Step 2: Check for presence of a CoA-like moiety.
    # We search for part of the adenosine moiety that is typical in CoA (e.g. adenine ring).
    # This SMARTS looks for an aromatic bicyclic ring containing nitrogen.
    coa_smarts = "n1cnc2"  
    coa_pat = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(coa_pat):
        return False, "No CoA moiety detected (adenine fragment missing)"
        
    # Step 3: Isolate the fatty acyl chain.
    # The fatty acyl chain is attached to the carbonyl carbon. We perform a DFS starting at the carbonyl carbon
    # and travel only through carbon atoms. We do not cross the bond to the sulfur or travel into the carbonyl oxygen.
    acyl_chain_atoms = set()
    visited = set()
    
    def dfs(atom_idx):
        visited.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        # Only consider carbon atoms.
        if atom.GetAtomicNum() == 6:
            acyl_chain_atoms.add(atom_idx)
        # Traverse neighbors, except we do not go from the carbonyl carbon into the sulfur or the carbonyl oxygen.
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            # Skip if already visited.
            if nbr_idx in visited:
                continue
            # At the starting carbonyl, do not cross into the S or the double-bonded O.
            if atom_idx == carbonyl_idx:
                # Get bond properties to decide direction.
                bond = mol.GetBondBetweenAtoms(atom_idx, nbr_idx)
                # If the neighbor is the sulfur, skip it.
                if nbr_idx == sulfur_idx:
                    continue
                # Skip any oxygen that is bonded by a double bond (the carbonyl oxygen).
                if nbr.GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() == 2.0:
                    continue
            # For subsequent atoms, only traverse if the neighbor is carbon.
            if nbr.GetAtomicNum() == 6:
                dfs(nbr_idx)
    
    # Begin DFS from the carbonyl carbon.
    dfs(carbonyl_idx)
    
    # The chain length is the number of carbon atoms found in this connected acyl fragment.
    chain_length = len(acyl_chain_atoms)
    
    if chain_length <= 22:
        return False, f"Fatty acyl chain has {chain_length} carbons; must be greater than 22"
    else:
        return True, f"Fatty acyl chain has {chain_length} carbons, which is >22"
        
# (Optional) To test the function, you can call it with one of the provided SMILES strings.
# For example:
# result, reason = is_very_long_chain_fatty_acyl_CoA("CC\\C=C/C\\C=C/C(=O)SCCNC(=O)CCNC(=O)[C@H](O)...")
# print(result, reason)