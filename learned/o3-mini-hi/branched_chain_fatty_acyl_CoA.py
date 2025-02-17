"""
Classifies: CHEBI:61912 branched-chain fatty acyl-CoA
"""
"""
Classifies: branched-chain fatty acyl-CoA

A branched-chain fatty acyl-CoA is defined as a fatty acyl-CoA derived from a branched-chain fatty acid 
(i.e. a fatty acyl group that contains at least one branch in its carbon chain) linked through a thioester bond 
(i.e. –C(=O)S–) to the Coenzyme A moiety.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_branched_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acyl-CoA based on its SMILES string.
    The decision is based on two steps:
      1. Checking for the presence of a CoA moiety using a canonical substructure pattern.
      2. Locating a thioester group ([#6](=O)S) and validating that the acyl (fatty acid) part 
         contains at least one branch (i.e. at least one carbon atom in the acyl fragment is connected 
         to three or more carbon atoms).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a branched-chain fatty acyl-CoA, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Step 1: Check for the Coenzyme A (CoA) moiety.
    # We look for a common fragment found in CoA esters.
    # Here we use a fragment pattern: SCCNC(=O)CCNC(=O)
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not found"

    # Step 2: Find the thioester group: a carbonyl carbon attached to S.
    # SMARTS: a carbon (atomic number 6) with a double-bonded O and a single bond to S.
    thioester_pattern = Chem.MolFromSmarts("[#6](=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester group (acyl-CoA linkage) found"

    # Helper function: given a starting carbon, traverse the fatty acyl chain (only through carbon atoms)
    # while not going back through the carbonyl (parent) or into CoA.
    def get_chain_atoms(start_idx, exclude_idx):
        """
        Traverse the chain from the starting atom (start_idx). Do not traverse into the atom given in exclude_idx.
        Only follow bonds between carbon atoms.
        Returns a set of atom indices for the acyl chain.
        """
        chain_atoms = set()
        frontier = [start_idx]
        while frontier:
            current = frontier.pop()
            if current in chain_atoms:
                continue
            chain_atoms.add(current)
            atom = mol.GetAtomWithIdx(current)
            # Only follow carbon neighbors (atomic number 6)
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                # Do not go back to the excluded atom (usually the carbonyl carbon)
                if nbr_idx == exclude_idx:
                    continue
                if nbr.GetAtomicNum() == 6 and nbr_idx not in chain_atoms:
                    frontier.append(nbr_idx)
        return chain_atoms

    # Now check each thioester match.
    # In the match tuple for "[#6](=O)S", the order is:
    # match[0] = carbonyl carbon, match[1] = oxygen (from =O), match[2] = sulfur (from S).
    branched_found = False
    for match in thioester_matches:
        carbonyl_idx = match[0]
        sulfur_idx = match[2]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        
        # Select the neighbor of the carbonyl carbon that is NOT the oxygen (already in match)
        acyl_root = None
        for nbr in carbonyl_atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            # Skip the sulfur (the thioester link) and the oxygen (double-bonded) 
            if nbr_idx in (match[1], sulfur_idx):
                continue
            # We assume the acyl chain is attached via a carbon atom.
            if nbr.GetAtomicNum() == 6:
                acyl_root = nbr_idx
                break
        if acyl_root is None:
            continue  # No valid acyl chain found for this thioester; try next
        
        # Collect all carbon atoms belonging to the fatty acyl chain.
        acyl_chain = get_chain_atoms(acyl_root, carbonyl_idx)
        # For a fatty acid chain, we would expect at least a few carbons.
        if len(acyl_chain) < 3:
            continue
        
        # Check for a branching point within the acyl chain:
        # For each atom in the chain, count its neighbors that are in the chain.
        # In a linear chain, internal carbons typically have 2 neighbors and terminal carbons 1.
        # A branching point will have degree 3 or more within the chain.
        for atom_idx in acyl_chain:
            atom = mol.GetAtomWithIdx(atom_idx)
            neighbor_count = sum(1 for nbr in atom.GetNeighbors() if nbr.GetIdx() in acyl_chain)
            if neighbor_count >= 3:
                branched_found = True
                break
        if branched_found:
            break

    if not branched_found:
        return False, "No branched fatty acyl chain detected"

    return True, "Contains a thioester-linked fatty acyl chain with a branch and a CoA moiety"

# Example usage:
if __name__ == "__main__":
    # Test with one of the provided examples: 2-methylbutanoyl-CoA
    example_smiles = "CCC(C)C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
    result, reason = is_branched_chain_fatty_acyl_CoA(example_smiles)
    print(result, reason)