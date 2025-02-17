"""
Classifies: CHEBI:25676 oligopeptide
"""
"""
Classifies: Oligopeptide – a peptide containing a relatively small number of amino acids.

This improved heuristic uses several steps:
  - It searches for peptide (amide) bonds using the SMARTS "C(=O)N".
  - It estimates the number of amino acid residues as (number of peptide bonds + 1) and requires that
    count to be between 2 and 10.
  - It counts occurrences of chiral alpha–carbon units (using [C@H](N) and [C@@H](N)).
  - It then verifies that all peptide–bond atoms (the carbonyl carbon and the adjacent amide nitrogen)
    are in a single contiguous fragment and further that they form a continuous linear chain.
  - To do this, for each peptide bond we add an undirected edge between the carbonyl carbon and amide nitrogen.
    In a correct peptide (with n bonds) the longest simple path in this subgraph should have exactly n edges.
  - It checks that the overall molecular weight roughly matches the residue count
    (using 60–170 Da per residue).

Note: This heuristic has limitations and may miss some peptides or include borderline molecules.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_oligopeptide(smiles: str):
    """
    Determines whether the given molecule (via its SMILES) is a small peptide (oligopeptide).
    
    Heuristics used:
      - Search for peptide (amide) bonds using SMARTS "C(=O)N".
      - Estimate residue count = (# peptide bonds + 1) and require 2 to 10 residues.
      - Count alpha–carbon centers using SMARTS "[C@H](N)" and "[C@@H](N)".
         * For dipeptides, require exactly 2 alpha centers.
         * For longer peptides, require at least 2.
      - Verify that all peptide bond atoms are in a single fragment.
      - Build an undirected graph from peptide bond atoms (carbonyl carbon and adjacent nitrogen)
        and require that the longest simple linear path (by edges) is equal to the number of peptide bonds.
      - Check that the molecular weight falls roughly in the range expected (60–170 Da per residue).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an oligopeptide, False otherwise.
        str: An explanation of the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find peptide (amide) bonds with SMARTS pattern "C(=O)N".
    peptide_bond_smarts = Chem.MolFromSmarts("C(=O)N")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_smarts)
    n_peptide_bonds = len(peptide_bond_matches)
    
    if n_peptide_bonds == 0:
        return False, "No peptide (amide) bonds found"
    
    # Estimate the number of residues.
    n_residues = n_peptide_bonds + 1
    if n_residues < 2:
        return False, "Too few peptide bonds to form a peptide (need at least 2 residues)"
    if n_residues > 10:
        return False, f"Found {n_residues} amino acid residues which is too many for an oligopeptide"
    
    # Count occurrences of alpha–carbon atoms (with chiral tag) using SMARTS patterns.
    alpha_smarts1 = Chem.MolFromSmarts("[C@H](N)")
    alpha_smarts2 = Chem.MolFromSmarts("[C@@H](N)")
    alpha_matches1 = mol.GetSubstructMatches(alpha_smarts1)
    alpha_matches2 = mol.GetSubstructMatches(alpha_smarts2)
    alpha_atoms = set()
    for match in alpha_matches1:
        # match[0] is the chiral carbon.
        alpha_atoms.add(match[0])
    for match in alpha_matches2:
        alpha_atoms.add(match[0])
    n_alpha = len(alpha_atoms)
    
    # For a dipeptide, require exactly 2 alpha–carbon centers.
    if n_residues == 2 and n_alpha != 2:
        return False, f"For a dipeptide, exactly 2 alpha–carbon centers are expected; found {n_alpha}"
    if n_residues >= 3 and n_alpha < 2:
        return False, f"Expected at least 2 alpha–carbon centers; found {n_alpha}"
    
    # Ensure that all atoms involved in the peptide bonds (the carbonyl carbon and the adjacent nitrogen)
    # are in a single contiguous fragment.
    pb_atoms = set()
    # Each match (from SMARTS "C(=O)N") returns three atoms: index0 = C(carbonyl), index1 = O, index2 = adjacent N.
    for match in peptide_bond_matches:
        pb_atoms.add(match[0])  # carbonyl carbon
        pb_atoms.add(match[2])  # amide nitrogen
    
    fragments = Chem.GetMolFrags(mol, asMols=False)
    n_frag_with_pb = sum(1 for frag in fragments if any(atom in frag for atom in pb_atoms))
    if n_frag_with_pb > 1:
        return False, "Peptide bond atoms are not contained in a single contiguous fragment"
    
    # Build an undirected graph from peptide bond atoms.
    # For each peptide bond match, add an edge between the carbonyl carbon and the amide nitrogen.
    graph = {}
    def add_edge(a, b):
        graph.setdefault(a, set()).add(b)
        graph.setdefault(b, set()).add(a)
    
    for match in peptide_bond_matches:
        c_atom = match[0]
        n_atom = match[2]
        add_edge(c_atom, n_atom)
    
    # Check that the peptide bonds form a continuous linear chain:
    # In a truly linear chain with n peptide bonds there will be a simple path of length n (in edges).
    # We use DFS to find the longest simple path within the graph.
    longest_path = 0
    def dfs(node, visited, length):
        nonlocal longest_path
        longest_path = max(longest_path, length)
        for neighbor in graph.get(node, []):
            if neighbor not in visited:
                dfs(neighbor, visited | {neighbor}, length + 1)
                
    for node in graph:
        dfs(node, {node}, 0)
    if longest_path < n_peptide_bonds:
        return False, (f"Peptide bond connectivity is not linear; "
                       f"expected a chain of {n_peptide_bonds} bonds but longest path found has {longest_path} bonds")
    
    # Calculate additional properties.
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Roughly expected molecular weight bounds.
    min_expected = n_residues * 60
    max_expected = n_residues * 170
    if mol_wt < min_expected:
        return False, (f"Molecular weight ({mol_wt:.1f} Da) is too low for a {n_residues}-residue peptide "
                       f"(expected at least {min_expected} Da)")
    if mol_wt > max_expected:
        return False, (f"Molecular weight ({mol_wt:.1f} Da) is too high for a {n_residues}-residue peptide "
                       f"(expected at most {max_expected} Da)")
    
    reason = (f"Detected {n_peptide_bonds} peptide bond(s) (≈{n_residues} residue(s)), {n_rotatable} rotatable bond(s), "
              f"{n_alpha} alpha–carbon center(s), and MW of {mol_wt:.1f} Da. The peptide-bond atoms form a "
              "single contiguous linear chain. This is consistent with an oligopeptide.")
    
    return True, reason

# For local testing one might do:
if __name__ == "__main__":
    # Example: Leu-Trp dipeptide
    test_smiles = "CC(C)C[C@H](N)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(O)=O"
    result, expl = is_oligopeptide(test_smiles)
    print(result, expl)