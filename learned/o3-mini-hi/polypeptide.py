"""
Classifies: CHEBI:15841 polypeptide
"""
#!/usr/bin/env python
"""
Classifies: Polypeptide – defined as a peptide containing ten or more amino acid residues.
This program detects candidate peptide bonds by checking for a carbonyl (C double‐bonded to O) 
bonded to a nitrogen that has an extra carbon neighbor. It then builds a graph connecting candidate
peptide bonds and performs a DFS search with cycle detection (i.e. tracking visited nodes in the
current path) to determine the longest continuous chain.
For linear peptides (both free N–terminus and free C–terminus), the residue count is taken as (# peptide bonds + 1);
for cyclic or blocked peptides, the residue count is approximated as the number of peptide bonds.
A molecule is classified as a polypeptide if it has at least 10 residues.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polypeptide(smiles: str):
    """
    Determines if a molecule is a polypeptide (>=10 amino acid residues) based on its SMILES string.
    
    The algorithm:
      1. Parse the SMILES string.
      2. Identify candidate peptide bonds. For each bond between carbon and nitrogen,
         the carbon must have a double bond to at least one oxygen (a carbonyl),
         and the nitrogen must be bonded to at least one carbon other than the one in the bond.
      3. Record candidate bonds as tuples: (carbon_idx, nitrogen_idx).
      4. Build a directed graph linking these candidate bonds. For each candidate bond (C, N),
         search the neighbors of the nitrogen (excluding the paired carbon) for a possible alpha carbon.
         If that alpha carbon is bonded to the carbon atom of another candidate peptide bond,
         add an edge from the first candidate to that candidate.
      5. Traverse the graph using DFS with cycle detection (to avoid infinite recursion)
         to compute the longest simple chain of candidate peptide bonds.
      6. Examine the molecule for free terminal groups (free amine and free acid).
         For a linear peptide chain, the residue count is (# peptide bonds + 1); otherwise, it is the number of bonds.
      7. If the estimated residue count is at least 10, return True; otherwise, return False.
    
    Args:
       smiles (str): SMILES string of the molecule.
    
    Returns:
       (bool, str): (True, reason) if the molecule appears to be a polypeptide, else (False, reason).
    """
    # Step 1: Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 2: Identify candidate peptide bonds.
    peptide_bonds = []  # Each entry is a tuple (c_idx, n_idx)
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        # Check if one atom is carbon and the other is nitrogen.
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 7:
            c_atom = a1
            n_atom = a2
        elif a1.GetAtomicNum() == 7 and a2.GetAtomicNum() == 6:
            n_atom = a1
            c_atom = a2
        else:
            continue

        # Require that the carbon has at least one double bond to oxygen (carbonyl).
        has_carbonyl = False
        for nbr in c_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                bond_CO = mol.GetBondBetweenAtoms(c_atom.GetIdx(), nbr.GetIdx())
                if bond_CO is not None and bond_CO.GetBondTypeAsDouble() >= 2.0:
                    has_carbonyl = True
                    break
        if not has_carbonyl:
            continue

        # Require that the nitrogen is connected to at least one additional carbon (besides c_atom).
        has_extra_c = False
        for nbr in n_atom.GetNeighbors():
            if nbr.GetIdx() == c_atom.GetIdx():
                continue
            if nbr.GetAtomicNum() == 6:
                has_extra_c = True
                break
        if not has_extra_c:
            continue

        peptide_bonds.append((c_atom.GetIdx(), n_atom.GetIdx()))
    
    if not peptide_bonds:
        return False, "No candidate peptide bonds found"

    # Step 4: Build the peptide bond graph.
    # We'll represent each candidate peptide bond by its index in the 'peptide_bonds' list.
    bond_graph = { i: [] for i in range(len(peptide_bonds)) }
    # Map from carbon index to a list of candidate bond indices where that atom is the carbon involved.
    c_to_bond = {}
    for i, (c_idx, n_idx) in enumerate(peptide_bonds):
        c_to_bond.setdefault(c_idx, []).append(i)
    
    # For each candidate bond, try to link to a following candidate bond.
    for i, (prev_c_idx, n_idx) in enumerate(peptide_bonds):
        n_atom = mol.GetAtomWithIdx(n_idx)
        # Look at neighbors of the nitrogen (excluding the paired carbon).
        for nbr in n_atom.GetNeighbors():
            if nbr.GetIdx() == prev_c_idx:
                continue
            if nbr.GetAtomicNum() != 6:
                continue
            alpha_idx = nbr.GetIdx()
            # Check if this "alpha carbon" is connected to the carbon of another peptide bond.
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
    
    # Step 5: Find the longest chain using DFS with cycle detection.
    # We define a DFS function that passes a set 'visited' (of candidate indices already in the current path).
    def dfs(node, visited):
        # Avoid cycles: if node is already in the current path, return 0.
        if node in visited:
            return 0
        max_length = 1  # At least the current candidate bond
        visited.add(node)
        for nxt in bond_graph.get(node, []):
            chain_length = 1 + dfs(nxt, visited)
            if chain_length > max_length:
                max_length = chain_length
        visited.remove(node)
        return max_length

    longest_chain = 0
    for i in range(len(peptide_bonds)):
        chain_length = dfs(i, set())
        if chain_length > longest_chain:
            longest_chain = chain_length

    # Step 6: Check for free terminal groups to classify the peptide as linear or not.
    free_amine = Chem.MolFromSmarts("[NH2]")
    free_acid = Chem.MolFromSmarts("C(=O)[OH]")
    has_free_amine = mol.HasSubstructMatch(free_amine)
    has_free_acid = mol.HasSubstructMatch(free_acid)
    
    if has_free_amine and has_free_acid:
        peptide_type = "linear"
        estimated_residues = longest_chain + 1  # For linear peptides, residue count = peptide bonds + 1.
    else:
        peptide_type = "cyclic or blocked"
        estimated_residues = longest_chain  # Approximation for cyclic/blocked peptides.
    
    if estimated_residues >= 10:
        return True, f"Detected {estimated_residues} amino acid residues (peptide appears {peptide_type})."
    else:
        return False, f"Detected only {estimated_residues} amino acid residues (need at least 10)."

# Example usage:
if __name__ == "__main__":
    # A few example SMILES strings (one from the supplied list is given here for testing).
    examples = {
        "nisin": "CC[C@H](C)[C@H](N)C(=O)N\\C(=C/C)C(=O)N[C@@H]1CSC[C@H](NC(=O)[C@H](CC(C)C)NC(=O)C(=C)NC(=O)[C@@H](NC1=O)[C@@H](C)CC)C(=O)N[C@@H]1[C@H](C)SC[C@H](NC(=O)CNC(=O)[C@@H]2CCCN2C1=O)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H]1[C@H](C)SC[C@H](NC(=O)CNC(=O)[C@H](CCSC)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C)NC(=O)CNC1=O)C(=O)N[C@@H](CC(N)=O)C(=O)N[C@@H](CCSC)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)C(N)=O"
    }
    
    for name, sm in examples.items():
        status, reason = is_polypeptide(sm)
        print("NAME:", name)
        print("Is polypeptide?", status)
        print("Reason:", reason)
        print()