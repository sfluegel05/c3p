"""
Classifies: CHEBI:15841 polypeptide
"""
Below is one possible improved solution. In our reasoning we noticed that simply “counting” any C–N bond that meets a minimal carbonyl–NH condition (as done before) may count bonds from fragments that are not connected in a single peptide‐backbone. Instead, we now try to “chain‐link” peptide bonds into a continuous backbone. In our improved approach we first flag bonds that look like peptide bonds (i.e. a C atom that is double‐bonded to oxygen and connected to a N atom that itself is attached to at least one other carbon) and record each as a (carbon,index, nitrogen_index) tuple. Then we build a (directed) “backbone graph” in which one peptide bond “leads” to another if the nitrogen of the first (after its “alpha carbon” connection) is attached to the carbonyl carbon of the next. (In a linear peptide the “alpha carbon” is the one attached to the amide N atom aside from the C=O.) Finally we search for the longest path (i.e. chain of linked peptide bonds) and then adjust the residue count depending on whether free termini (a free primary amine and carboxylic acid) are present. (For a linear chain the estimated number of amino acid residues equals the peptide‐bond count + 1; otherwise we assume a cyclic or blocked chain where the residue count roughly equals the number of peptide bonds.) This “chain‐detection” heuristic better discriminates molecules that truly have a continuous peptide backbone from false positives that simply contain one or two amide bonds only. (Of course, the heuristic is still imperfect.) 

Below is the complete Python program:

------------------------------------------------

#!/usr/bin/env python
"""
Classifies: Polypeptide – defined as a peptide containing ten or more amino acid residues.
This version attempts to improve on the earlier heuristic by (a) first detecting bonds that appear to be peptide bonds 
(using a carbonyl–nitrogen rule) and then (b) linking them into a continuous peptide backbone.
If a free N‐terminus ([NH2]) and free C‐terminus (C(=O)[OH]) are found we assume a linear peptide.
Otherwise we assume a cyclic or blocked structure so that the peptide “residue count” equals the number of peptide bonds.
Note: This is a heuristic still and may be fooled in some edge‐cases.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polypeptide(smiles: str):
    """
    Determines if a molecule is a polypeptide (≥10 amino acid residues) based on its SMILES string.
    The algorithm:
      1. Parse the SMILES.
      2. Loop over bonds and mark each bond between a carbon and a nitrogen as a candidate peptide bond
         if:
           (a) the carbon has at least one neighbor oxygen with which it forms a double bond (C=O)
           (b) the nitrogen is connected to at least one other carbon (aside from the bonded carbon).
      3. For each candidate peptide bond, record the (carbon_index, nitrogen_index).
      4. Build a directed graph over these peptide bonds to try to link bonds into a continuous backbone.
         For each peptide bond bond_i = (C_i, N_i), we look at the neighbors of its nitrogen (other than C_i)
         that are carbon atoms. If one such neighbor (call it CA) is itself the carbon of a candidate peptide bond bond_j,
         then we say bond_i “leads” to bond_j. (The idea is that in an amino acid residue the amide nitrogen is attached to an alpha carbon that, in turn, is attached to the carbonyl of the next residue.)
      5. Find the longest chain (i.e. the longest directed path) in this graph.
         Then, if a free amine and free acid can be detected (using SMARTS patterns) we consider the chain linear and assign
         residue_count = (longest chain peptide bond count + 1); otherwise, we use residue_count = longest chain peptide bond count.
      6. Return True if residue_count ≥ 10.
    
    Args:
       smiles (str): SMILES string of the molecule.
    
    Returns:
       (bool, str): (True, reason) if the molecule appears to be a polypeptide, else (False, reason).
    """
    # Parse molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Step 2: Identify candidate peptide bonds.
    # Each bond must be between a carbon and a nitrogen.
    # Condition: the carbon must be double-bonded to at least one oxygen,
    # and the nitrogen must have a neighboring carbon other than the one in the bond.
    peptide_bonds = []  # list of tuples: (c_idx, n_idx)
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        # Identify candidate C and N
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 7:
            c_atom = a1
            n_atom = a2
        elif a1.GetAtomicNum() == 7 and a2.GetAtomicNum() == 6:
            n_atom = a1
            c_atom = a2
        else:
            continue

        # Check carbon: must have at least one double bond to oxygen (carbonyl)
        has_carbonyl = False
        for nbr in c_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                b = mol.GetBondBetweenAtoms(c_atom.GetIdx(), nbr.GetIdx())
                # RDKit represents bond order as a float via GetBondTypeAsDouble()
                if b is not None and b.GetBondTypeAsDouble() >= 2.0:
                    has_carbonyl = True
                    break
        if not has_carbonyl:
            continue

        # Check nitrogen: must be connected to at least one other carbon (besides c_atom)
        has_other_c = False
        for nbr in n_atom.GetNeighbors():
            if nbr.GetIdx() == c_atom.GetIdx():
                continue
            if nbr.GetAtomicNum() == 6:
                has_other_c = True
                break
        if not has_other_c:
            continue

        peptide_bonds.append( (c_atom.GetIdx(), n_atom.GetIdx()) )

    # If no peptide bonds were found, we cannot be a polypeptide.
    if not peptide_bonds:
        return False, "No candidate peptide bonds found"

    # Step 4: Build a graph of peptide bonds.
    # We want to link bond i -> bond j if:
    #   The nitrogen of bond i (n_i) connects to an alpha carbon that is the carbon (c_j) of bond j.
    # Here we keep an index for each candidate peptide bond.
    bond_graph = { i: [] for i in range(len(peptide_bonds)) }
    # Build a mapping from carbon atom index to peptide bond indices for which it is the carbon partner.
    c_to_bond = {}
    for i, (c_idx, n_idx) in enumerate(peptide_bonds):
        c_to_bond.setdefault(c_idx, []).append(i)
    # Now, for each peptide bond i, inspect the nitrogen side.
    for i, (c_idx, n_idx) in enumerate(peptide_bonds):
        n_atom = mol.GetAtomWithIdx(n_idx)
        # Look for an adjacent carbon (potential alpha carbon) from the nitrogen (but not back to c_idx)
        for nbr in n_atom.GetNeighbors():
            if nbr.GetIdx() == c_idx:
                continue
            if nbr.GetAtomicNum() != 6: 
                continue
            alpha_idx = nbr.GetIdx()
            # Now, if this alpha carbon is itself directly bonded to a carbon that is the carbon component of another peptide bond, add an edge.
            for nbr2 in nbr.GetNeighbors():
                if nbr2.GetIdx() == n_idx:
                    continue
                if nbr2.GetAtomicNum() == 6:
                    # Check if nbr2 is in our mapping from carbon to peptide bond.
                    if nbr2.GetIdx() in c_to_bond:
                        for j in c_to_bond[nbr2.GetIdx()]:
                            # Avoid self-loop.
                            if j == i:
                                continue
                            bond_graph[i].append(j)
    # Step 5: Find the longest chain.
    # We use a DFS with memoization on our bond_graph nodes.
    memo = {}
    def dfs(node):
        if node in memo:
            return memo[node]
        max_length = 1  # counting the current peptide bond as length 1
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

    # Step 6: Determine if the chain appears linear.
    free_amine = Chem.MolFromSmarts("[NH2]")
    free_acid = Chem.MolFromSmarts("C(=O)[OH]")
    has_free_amine = mol.HasSubstructMatch(free_amine)
    has_free_acid = mol.HasSubstructMatch(free_acid)
    if has_free_amine and has_free_acid:
        peptide_type = "linear"
        estimated_residues = longest_chain + 1  # linear: number of bonds +1
    else:
        peptide_type = "cyclic or with blocked termini"
        estimated_residues = longest_chain  # cyclic: residue count roughly equals peptide bond count

    # Final classification: require ≥ 10 amino acid residues.
    if estimated_residues >= 10:
        return True, f"Detected {estimated_residues} amino acid residues (peptide appears {peptide_type})."
    else:
        return False, f"Detected only {estimated_residues} amino acid residues (need at least 10)."

# Example usage:
if __name__ == "__main__":
    # Example molecules: (only two examples shown; more may be tested)
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
        
"""
Reasoning summary:
Our earlier version simply tallied peptide bonds meeting minimal criteria. However, many false positives/negatives occurred 
because the bonds were not linked into a single continuous chain. Here we build a graph that “chains” candidate peptide bonds 
via the connecting (alpha) carbon from the amide nitrogen. In addition, we check for free terminal groups so that linear peptides 
count as bond_count+1 residues. Although still heuristic (and noting that extreme cases – cyclic peptides, branched or modified 
backbones – may be mis‐interpreted), this approach improves the detection of polypeptides with ten or more amino acid residues.
"""
  

