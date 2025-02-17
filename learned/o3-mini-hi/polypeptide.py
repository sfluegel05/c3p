"""
Classifies: CHEBI:15841 polypeptide
"""
"""
Classifies: Polypeptide – defined as a peptide containing ten or more amino acid residues.
This improved version tightens candidate peptide bond detection by requiring that the carbon
involved is sp2 (to support a carbonyl) and that the nitrogen is linked to a chiral (alpha) carbon.
It then builds a backbone graph from these candidate bonds by linking via the alpha carbon.
For linear peptides (with free N- and C-termini) the residue count is (# peptide bonds + 1);
for cyclic or blocked peptides, we approximate it as the number of bonds.
A molecule is classified as a polypeptide if it has at least 10 residues.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdchem import HybridizationType, ChiralType

def is_polypeptide(smiles: str):
    """
    Determines if a molecule is a polypeptide (>=10 amino acid residues) based on its SMILES string.
    
    The algorithm:
      1. Parse the SMILES.
      2. Identify candidate amide bonds. For a bond between a carbon and a nitrogen:
         - The carbon must be sp2-hybridized and have a double bond to an oxygen (i.e. a carbonyl).
         - The nitrogen must have at least one neighboring carbon that is chiral (candidate for an alpha carbon).
      3. Record each candidate as (c_idx, n_idx, alpha_idx) where alpha_idx is the neighbor chiral carbon.
         (If more than one chiral carbon is found, we take the first.)
      4. Build a directed graph where an edge goes from candidate bond i to candidate bond j
         if the alpha carbon from bond i is directly bonded to the carbon (c_idx) of candidate bond j.
      5. Use DFS (with cycle–avoidance) to compute the longest simple chain.
      6. Check the molecule for free terminal groups. For a linear peptide chain,
         the residue count is (# peptide bonds + 1); otherwise, it is the number of bonds.
      7. If the estimated residue count is at least 10, return True; else return False.
    
    Args:
       smiles (str): SMILES string of the molecule.
    
    Returns:
       tuple(bool, str): (True, reason) if the molecule appears to be a polypeptide,
                         else (False, reason).
    """
    # Step 1: Parse the SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 2: Identify candidate peptide bonds.
    # Each candidate will be a tuple: (c_idx, n_idx, alpha_idx)
    # where c_idx = index of the carbonyl carbon,
    #       n_idx = index of amide nitrogen,
    #       alpha_idx = index of the α–carbon (a chiral carbon attached to the N, if any)
    candidates = []
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if (a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 7):
            c_atom = a1
            n_atom = a2
        elif (a1.GetAtomicNum() == 7 and a2.GetAtomicNum() == 6):
            n_atom = a1
            c_atom = a2
        else:
            continue

        # Check that the carbon (c_atom) is sp2 (as expected for a carbonyl)
        if c_atom.GetHybridization() != HybridizationType.SP2:
            continue

        # Check that c_atom has at least one double bond to oxygen.
        has_carbonyl = False
        for nbr in c_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                bond_CO = mol.GetBondBetweenAtoms(c_atom.GetIdx(), nbr.GetIdx())
                # Check for a double bond (in RDKit bond types, GetBondTypeAsDouble() returns 2.0 for a double bond)
                if bond_CO and bond_CO.GetBondTypeAsDouble() >= 2.0:
                    has_carbonyl = True
                    break
        if not has_carbonyl:
            continue

        # Check that the nitrogen (n_atom) has at least one neighboring carbon that is chiral 
        # (candidate for the alpha carbon). Exclude the candidate carbon we already used.
        alpha_idx = None
        for nbr in n_atom.GetNeighbors():
            if nbr.GetIdx() == c_atom.GetIdx():
                continue
            if nbr.GetAtomicNum() == 6 and nbr.GetChiralTag() != ChiralType.CHI_UNSPECIFIED:
                alpha_idx = nbr.GetIdx()
                break
        if alpha_idx is None:
            continue

        candidates.append((c_atom.GetIdx(), n_atom.GetIdx(), alpha_idx))
    
    if not candidates:
        return False, "No candidate peptide bonds found"

    # Step 3: Build the peptide bond graph.
    # We assign each candidate bond an index. We want an edge from candidate i to candidate j
    # if the alpha carbon of candidate i is directly bonded to the carbon of candidate j.
    bond_graph = { i: [] for i in range(len(candidates)) }
    # Make a mapping from a carbon index (the carbonyl carbon) to candidate indices.
    c_to_candidate = {}
    for i, (c_idx, _, _) in enumerate(candidates):
        c_to_candidate.setdefault(c_idx, []).append(i)

    for i, (c_idx, n_idx, alpha_idx) in enumerate(candidates):
        # Get the alpha atom (from the candidate's N neighbors)
        alpha_atom = mol.GetAtomWithIdx(alpha_idx)
        # Look for neighbors of this alpha atom. If one of these is the carbon (c_idx) of another candidate,
        # then add an edge from candidate i -> candidate j.
        for nbr in alpha_atom.GetNeighbors():
            cand_c_idx = nbr.GetIdx()
            if cand_c_idx in c_to_candidate:
                for j in c_to_candidate[cand_c_idx]:
                    # Avoid self-linking.
                    if i == j:
                        continue
                    bond_graph[i].append(j)
    
    # Step 4: Find the longest chain using DFS with cycle detection.
    def dfs(node, visited):
        if node in visited:
            return 0  # cycle detected, do not count further
        max_length = 1  # count current candidate
        visited.add(node)
        for nxt in bond_graph.get(node, []):
            length_here = 1 + dfs(nxt, visited)
            if length_here > max_length:
                max_length = length_here
        visited.remove(node)
        return max_length

    longest_chain = 0
    for i in range(len(candidates)):
        chain_length = dfs(i, set())
        if chain_length > longest_chain:
            longest_chain = chain_length

    # Step 5: Check for free terminal groups to guess if the peptide is linear.
    # A linear peptide is likely if both a free amine and a free acid are present.
    free_amine = Chem.MolFromSmarts("[NH2]")
    free_acid = Chem.MolFromSmarts("C(=O)[OH]")
    has_free_amine = mol.HasSubstructMatch(free_amine)
    has_free_acid = mol.HasSubstructMatch(free_acid)

    if has_free_amine and has_free_acid:
        peptide_type = "linear"
        estimated_residues = longest_chain + 1  # chain length (# bonds) + one terminal residue
    else:
        peptide_type = "cyclic or blocked"
        estimated_residues = longest_chain  # for cyclic or blocked peptides, we use number of bonds

    if estimated_residues >= 10:
        return True, f"Detected {estimated_residues} amino acid residues (peptide appears {peptide_type})."
    else:
        return False, f"Detected only {estimated_residues} amino acid residues (need at least 10)."

# Example usage:
if __name__ == "__main__":
    # Try a few example SMILES.
    example_smiles = {
        "nisin": "CC[C@H](C)[C@H](N)C(=O)N\\C(=C/C)C(=O)N[C@@H]1CSC[C@H](NC(=O)[C@H](CC(C)C)NC(=O)C(=C)NC(=O)[C@@H](NC1=O)[C@@H](C)CC)C(=O)N[C@@H]1[C@H](C)SC[C@H](NC(=O)CNC(=O)[C@@H]2CCCN2C1=O)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H]1[C@H](C)SC[C@H](NC(=O)CNC(=O)[C@H](CCSC)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C)NC(=O)CNC1=O)C(=O)N[C@@H](CC(N)=O)C(=O)N[C@@H](CCSC)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)C(N)=O"
    }
    
    for name, sm in example_smiles.items():
        status, reason = is_polypeptide(sm)
        print("NAME:", name)
        print("Is polypeptide?", status)
        print("Reason:", reason)
        print()