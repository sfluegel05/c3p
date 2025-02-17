"""
Classifies: CHEBI:83139 long-chain fatty acyl-CoA(4-)
"""
#!/usr/bin/env python
"""
Classifies: long-chain fatty acyl-CoA(4-)
Definition: A fatty acyl-CoA(4-) arising from deprotonation of the phosphate and diphosphate OH groups
of any long-chain fatty acyl-CoA; major species at pH 7.3.
This implementation verifies that the molecule contains (a) an adenine nucleotide fragment,
(b) a pantetheine fragment, (c) a thioester group linking a fatty acyl chain and (d) an acyl chain 
with a longest linear carbon path of 12–22 carbons.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines whether the given SMILES string belongs to the class long-chain fatty acyl-CoA(4-).

    Process:
      1. Parse the SMILES.
      2. Ensure that an adenine nucleotide fragment is present (via several SMARTS patterns).
      3. Verify the presence of a pantetheine fragment.
      4. Look for a thioester group ([C;!R](=O)[S]) which links the fatty acyl chain.
      5. For each thioester hit, select the carbonyl carbon and, from its neighbor (the acyl group),
         perform a DFS to compute the longest carbon–chain (all atoms are carbon) starting from that atom.
         This helps to allow minor branching (as in iso–acyl species) while still counting the main chain.
      6. Accept if at least one thioester leads to a longest chain path between 12 and 22 carbons.

    Args:
        smiles (str): SMILES of the molecule.

    Returns:
        (bool, str): A tuple with True and an explanation if the molecule fits the class,
                     otherwise False with an explanation.
    """
    # --- Step 1: Parse SMILES ---
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Step 2: Check for adenine fragment using alternative SMARTS patterns ---
    adenine_smarts = [
        "n1cnc2ncnc12",         # plain adenine core
        "n1cnc2c(n)ncnc12",      # variant with a methyl substituent
        "n1cnc2c(N)ncnc12"       # variant with an amino substituent
    ]
    found_adenine = False
    for smarts in adenine_smarts:
        patt = Chem.MolFromSmarts(smarts)
        if patt and mol.HasSubstructMatch(patt):
            found_adenine = True
            break
    if not found_adenine:
        return False, "Adenine nucleotide fragment not found; not a CoA derivative"
    
    # --- Step 3: Check for the pantetheine fragment ---
    pantetheine_patt = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if not pantetheine_patt or not mol.HasSubstructMatch(pantetheine_patt):
        return False, "Pantetheine fragment not found; not a CoA derivative"
    
    # --- Step 4: Locate a thioester group ---
    thioester_patt = Chem.MolFromSmarts("[C;!R](=O)[S]")
    thioester_matches = mol.GetSubstructMatches(thioester_patt)
    if not thioester_matches:
        return False, "No thioester group (fatty acid linkage) found"
    
    # --- Step 5: DFS to compute the longest carbon chain from the fatty acyl start ---
    def dfs_longest_chain(atom, parent_idx, visited):
        """
        Performs depth-first search from the given atom over neighboring carbons to find
        the length of the longest path (in number of carbon atoms). 'parent_idx' is used
        to avoid immediately backtracking.
        """
        max_length = 1  # count this atom itself
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:  # only follow carbon atoms
                continue
            if nbr.GetIdx() == parent_idx:
                continue
            if nbr.GetIdx() in visited:
                continue
            new_visited = visited | {nbr.GetIdx()}
            length = 1 + dfs_longest_chain(nbr, atom.GetIdx(), new_visited)
            if length > max_length:
                max_length = length
        return max_length

    MIN_CHAIN = 12
    MAX_CHAIN = 22
    acyl_chain_valid = False
    chain_length_found = 0

    # For each thioester match, try to assess the attached acyl chain.
    for match in thioester_matches:
        carbonyl_idx, sulfur_idx = match[0], match[1]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        fatty_start = None
        # From the carbonyl, pick a neighbor that is not the sulfur and is a carbon.
        for nbr in carbonyl_atom.GetNeighbors():
            if nbr.GetIdx() == sulfur_idx:
                continue
            if nbr.GetAtomicNum() == 6:
                fatty_start = nbr
                break
        if fatty_start is None:
            continue  # try next thioester
        
        # To avoid backtracking into the carbonyl, start DFS from fatty_start
        # Start with visited set containing fatty_start (its index) and treat carbonyl as parent.
        longest_path = 1 + dfs_longest_chain(fatty_start, carbonyl_atom.GetIdx(), {fatty_start.GetIdx()})
        if MIN_CHAIN <= longest_path <= MAX_CHAIN:
            acyl_chain_valid = True
            chain_length_found = longest_path
            break

    if not acyl_chain_valid:
        return False, f"Fatty acyl chain missing or longest linear carbon path is not within {MIN_CHAIN}-{MAX_CHAIN}."
    
    return True, (
        f"Contains CoA moiety (adenine and pantetheine fragments) and a thioester-linked fatty acyl chain "
        f"with longest linear carbon path of {chain_length_found} carbons, as expected for long-chain fatty acyl-CoA(4-)."
    )


# For testing purposes:
if __name__ == "__main__":
    # Example test (one true positive)
    test_smiles = (
        "CCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)"
        "[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)"
        "OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
    )
    result, explanation = is_long_chain_fatty_acyl_CoA_4__(test_smiles)
    print(result, explanation)