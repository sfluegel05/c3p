"""
Classifies: CHEBI:15841 polypeptide
"""
#!/usr/bin/env python
"""
Classifies: Polypeptide – defined as a peptide containing ten or more amino acid residues.
This program attempts to detect a continuous peptide backbone by first identifying candidate peptide bonds
(using a carbonyl-to-nitrogen rule) and then linking them into a backbone chain.
If both a free N-terminus ([NH2]) and a free C-terminus (C(=O)[OH]) are found, we treat the chain as linear
and count residues as (# peptide bonds + 1). Otherwise, for cyclic/blocked peptides we use the number of peptide bonds
as an estimate for the residue count.
Note: This approach is heuristic and may be imperfect for very complex/multicyclic molecules.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polypeptide(smiles: str):
    """
    Determines if a molecule is a polypeptide (>=10 amino acid residues) based on its SMILES string.
    The algorithm:
      1. Parse the SMILES string.
      2. Identify candidate peptide bonds. A bond between a carbon and a nitrogen is marked candidate if:
            - The carbon atom is double-bonded to at least one oxygen atom (i.e. a carbonyl).
            - The nitrogen atom is attached to at least one additional carbon (other than the one in the bond).
      3. Record candidate bonds as tuples: (carbon_idx, nitrogen_idx).
      4. Build a directed graph linking these bonds. For each candidate bond (C_i, N_i), we look at the neighbors
         of N_i (excluding C_i) for a carbon (the possible alpha carbon). If that neighbor is directly bonded to the carbon
         of another candidate peptide bond, we add a link from the first bond to that candidate.
      5. Traverse this graph (using DFS with memoization) to find the longest continuous chain.
      6. Check for free N-terminus and C-terminus using SMARTS. For linear peptides the estimated residue count is
         (# peptide bonds in longest chain + 1) while for cyclic/blocked peptides it equals the number of bonds.
      7. Return True if the estimated number of amino acid residues is at least 10.
    
    Args:
       smiles (str): SMILES string of the molecule.
    
    Returns:
       (bool, str): (True, reason) if the molecule appears to be a polypeptide, else (False, reason).
    """
    # Step 1: Parse molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 2: Identify candidate peptide bonds.
    peptide_bonds = []  # list of (c_idx, n_idx) pairs
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        # Identify if one atom is carbon and the other nitrogen.
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 7:
            c_atom = a1
            n_atom = a2
        elif a1.GetAtomicNum() == 7 and a2.GetAtomicNum() == 6:
            n_atom = a1
            c_atom = a2
        else:
            continue

        # Check carbon: must have a double bond to oxygen.
        has_carbonyl = False
        for nbr in c_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                bond_CO = mol.GetBondBetweenAtoms(c_atom.GetIdx(), nbr.GetIdx())
                if bond_CO is not None and bond_CO.GetBondTypeAsDouble() >= 2.0:
                    has_carbonyl = True
                    break
        if not has_carbonyl:
            continue

        # Check nitrogen: must be connected to at least one carbon other than c_atom.
        has_extra_c = False
        for nbr in n_atom.GetNeighbors():
            if nbr.GetIdx() == c_atom.GetIdx():
                continue
            if nbr.GetAtomicNum() == 6:
                has_extra_c = True
                break
        if not has_extra_c:
            continue

        peptide_bonds.append( (c_atom.GetIdx(), n_atom.GetIdx()) )
    
    if not peptide_bonds:
        return False, "No candidate peptide bonds found"
    
    # Step 4: Build the backbone graph.
    # Create a mapping from candidate index to list of candidate indices that might follow.
    bond_graph = { i: [] for i in range(len(peptide_bonds)) }
    # Map each carbon index to candidate bond indices where that carbon is the peptide bond carbon.
    c_to_bond = {}
    for i, (c_idx, n_idx) in enumerate(peptide_bonds):
        c_to_bond.setdefault(c_idx, []).append(i)
    
    # For each candidate peptide bond, look at the nitrogen's neighbors.
    for i, (prev_c_idx, n_idx) in enumerate(peptide_bonds):
        n_atom = mol.GetAtomWithIdx(n_idx)
        # Check neighbors of the nitrogen for a possible "alpha carbon" (exclude the carbon it is bonded to in this bond).
        for nbr in n_atom.GetNeighbors():
            if nbr.GetIdx() == prev_c_idx:
                continue
            if nbr.GetAtomicNum() != 6:
                continue
            alpha_idx = nbr.GetIdx()
            # Check if this alpha atom is next to another candidate peptide bond's carbon.
            for second_nbr in nbr.GetNeighbors():
                if second_nbr.GetIdx() == n_idx:
                    continue
                if second_nbr.GetAtomicNum() == 6:
                    cand_c_idx = second_nbr.GetIdx()
                    if cand_c_idx in c_to_bond:
                        for j in c_to_bond[cand_c_idx]:
                            if j == i:
                                continue
                            bond_graph[i].append(j)
    # Step 5: Find the longest chain using DFS with memoization.
    memo = {}
    def dfs(node):
        if node in memo:
            return memo[node]
        max_length = 1  # count current peptide bond
        for nxt in bond_graph.get(node, []):
            length = 1 + dfs(nxt)
            if length > max_length:
                max_length = length
        memo[node] = max_length
        return max_length

    longest_chain = 0
    for i in bond_graph:
        chain_length = dfs(i)
        if chain_length > longest_chain:
            longest_chain = chain_length

    # Step 6: Check for free terminal groups.
    free_amine = Chem.MolFromSmarts("[NH2]")
    free_acid = Chem.MolFromSmarts("C(=O)[OH]")
    has_free_amine = mol.HasSubstructMatch(free_amine)
    has_free_acid = mol.HasSubstructMatch(free_acid)
    
    if has_free_amine and has_free_acid:
        peptide_type = "linear"
        estimated_residues = longest_chain + 1  # linear peptides: bonds + 1
    else:
        peptide_type = "cyclic or blocked"
        estimated_residues = longest_chain  # approximate residue count equals bond count
    
    if estimated_residues >= 10:
        return True, f"Detected {estimated_residues} amino acid residues (peptide appears {peptide_type})."
    else:
        return False, f"Detected only {estimated_residues} amino acid residues (need at least 10)."

# Example usage:
if __name__ == "__main__":
    examples = {
        "nisin": "CC[C@H](C)[C@H](N)C(=O)N\\C(=C/C)C(=O)N[C@@H]1CSC[C@H](NC(=O)[C@H](CC(C)C)NC(=O)C(=C)NC(=O)[C@@H](NC1=O)[C@@H](C)CC)C(=O)N[C@@H]1[C@H](C)SC[C@H](NC(=O)CNC(=O)[C@@H]2CCCN2C1=O)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H]1[C@H](C)SC[C@H](NC(=O)CNC(=O)[C@H](CCSC)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C)NC(=O)CNC1=O)C(=O)N[C@@H](CC(N)=O)C(=O)N[C@@H](CCSC)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)C(N)=O",
        "Ac-Asp-N(6)-[2-(trifluoromethyl)benzoyl]-KATIGFEVQEE": "CC[C@H](C)[C@H](NC(=O)[C@@H](NC(=O)[C@H](C)NC(=O)[C@H](CCCCNC(=O)c1ccccc1C(F)(F)F)NC(=O)[C@H](CC(O)=O)NC(C)=O)[C@@H](C)O)C(=O)NCC(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](CCC(O)=O)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CCC(N)=O)C(=O)N[C@@H](CCC(O)=O)C(=O)N[C@@H](CCC(O)=O)C(O)=O"
    }
    
    for name, sm in examples.items():
        status, reason = is_polypeptide(sm)
        print("NAME:", name)
        print("Is polypeptide?", status)
        print("Reason:", reason)
        print()