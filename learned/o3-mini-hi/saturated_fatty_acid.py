"""
Classifies: CHEBI:26607 saturated fatty acid
"""
"""
Classifies: Saturated Fatty Acid
Definition: A fatty acid is “saturated” if it (i) contains a terminal carboxyl group,
(ii) has no carbon–carbon double or triple bonds outside of C=O groups,
and (iii) its alkyl (acyl) chain is “simple” – that is, aside from a possible terminal (ω)–OH and 
occasional carbonyl (C=O), any branching is limited to single methyl groups; notably, an α–carbon 
substituted by –OH or branches longer than methyl are not allowed.
Note: This heuristic‐based classifier was tuned to minimize false positives (e.g. complex diacids,
ether‐linked acids, or fatty acyl lipids) while still capturing “classical” (aliphatic) saturated fatty acids.
"""

from rdkit import Chem

def is_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a saturated fatty acid based on its SMILES string.
    
    Requirements (heuristic):
      1. The SMILES must be valid and contain only allowed elements (C, H, O, D).
      2. The molecule must be acyclic.
      3. There is exactly one carboxylic acid group (SMARTS "C(=O)[O;H,-]") 
         and its acid carbon is terminal (attached to exactly one carbon).
      4. The α–carbon (the single carbon attached to the carboxyl carbon) must not be substituted by oxygen.
      5. The molecule has no C=C or C≡C bonds between carbon atoms.
      6. The acyl “backbone” is extracted (via a DFS from the α–carbon through C–C bonds)
         and any branch (i.e. substituent off the main chain) must be just a methyl group.
      7. Any oxygen attached (by a single bond) to a backbone carbon is allowed only at the terminal ω–position.
         (Double–bonded carbonyl oxygens are permitted.)
         
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a saturated fatty acid, False otherwise.
        str: Explanation for the classification.
    """
    # Parse molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1: Check allowed elements: only C, H, O, and D (for deuterium)
    allowed = {"C", "H", "O", "D"}
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in allowed:
            return False, f"Disallowed element found ({atom.GetSymbol()}); not a simple fatty acid"
    
    # Step 2: Reject molecules with rings (to reduce complexity and avoid sugars or cyclic lipids)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Contains cyclic system; not a simple saturated fatty acid"
    
    # Step 3: Identify carboxylic acid group using SMARTS.
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H,-]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if len(acid_matches) != 1:
        return False, f"Expected exactly one carboxylic acid group, found {len(acid_matches)}"
    
    # Use the single match.
    match = acid_matches[0]
    acid_carbon_idx = match[0]  # acid carbon (the one double-bonded to O)
    acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
    
    # Step 4: Check that the acid carbon is terminal (has exactly one carbon neighbor)
    carbon_neighbors = [nbr for nbr in acid_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1:
        return False, "Carboxylic acid group is not terminal (acid carbon is linked to >1 carbon)"
    alpha = carbon_neighbors[0]  # α–carbon attached to acid moiety
    alpha_idx = alpha.GetIdx()
    
    # Rule: α–carbon must not be substituted by an -OH (or any oxygen by a single bond)
    for nbr in alpha.GetNeighbors():
        if nbr.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(alpha_idx, nbr.GetIdx()).GetBondType() == Chem.BondType.SINGLE:
            return False, "α–carbon is substituted by an –OH; not a typical saturated fatty acid"
    
    # Step 5: Check for any C=C or C≡C bonds (only allow C=O bonds)
    for bond in mol.GetBonds():
        bt = bond.GetBondType()
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if bt in [Chem.BondType.DOUBLE, Chem.BondType.TRIPLE]:
            # If both atoms are carbon then this is an unsaturation we do not allow.
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                return False, "Contains carbon–carbon unsaturation"
    
    # Create a set of indices for carbon atoms.
    carbon_ids = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6}
    
    # --- Helper: DFS to find the longest continuous carbon chain starting at 'start_idx'.
    def dfs_longest_chain(current_idx, visited):
        visited = visited | {current_idx}
        longest = [current_idx]
        for nbr in mol.GetAtomWithIdx(current_idx).GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            nbr_idx = nbr.GetIdx()
            if nbr_idx in visited:
                continue
            # Only traverse along carbon–carbon bonds.
            candidate = [current_idx] + dfs_longest_chain(nbr_idx, visited)
            if len(candidate) > len(longest):
                longest = candidate
        return longest
    
    # Start the chain from the α–carbon.
    main_chain = dfs_longest_chain(alpha_idx, set())
    if len(main_chain) < 1:
        return False, "Could not extract a carbon chain from the α–carbon"
    
    # For clarity, sort the main chain so that the acid side is excluded.
    # (In our DFS the order is not necessarily linear; we assume the DFS found one of the longest linear paths.)
    # We will treat main_chain as the “acyl chain” (not including the acid carbon).
    
    # --- Step 6: Check branching along the main chain.
    # The idea: for each carbon in the main chain, check for carbon neighbors that are not part of the linear chain.
    # They must be a terminal methyl (i.e. in the carbon graph, after removing the connection back, degree == 1).
    # Also, check oxygen substituents: allow them only if the bond is double (a carbonyl) or if the carbon is the terminal (ω–position).
    main_chain_set = set(main_chain)
    chain_len = len(main_chain)
    for i, atom_idx in enumerate(main_chain):
        atom = mol.GetAtomWithIdx(atom_idx)
        # Determine adjacent main chain indices (previous and next if they exist)
        adjacent_chain = set()
        if i > 0:
            adjacent_chain.add(main_chain[i-1])
        if i < chain_len - 1:
            adjacent_chain.add(main_chain[i+1])
        # Check each neighbor of this backbone carbon.
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            # Skip if neighbor is the acid carbon (only relevant at the chain terminus) or in the main chain.
            if nbr_idx == acid_carbon_idx or nbr_idx in adjacent_chain:
                continue
            # If neighbor is carbon (a branch)
            if nbr.GetAtomicNum() == 6:
                # In the subgraph of carbons, count how many connections nbr has excluding its bond back to 'atom'
                nbr_neighbors = [n.GetIdx() for n in nbr.GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIdx() != atom_idx]
                if len(nbr_neighbors) != 0:
                    # A branch longer than a methyl group.
                    return False, "Found branch substituent (longer than a methyl) on the acyl chain"
            # If neighbor is oxygen, check the bond type.
            elif nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(atom_idx, nbr_idx)
                # Allow if the oxygen is double-bonded (i.e., a carbonyl) OR if the carbon is the terminal ω–position.
                if bond.GetBondType() != Chem.BondType.DOUBLE:
                    # Allow only if this atom is the ω–end (last carbon in the main chain).
                    if i != chain_len - 1:
                        return False, "Backbone carbon (non–terminal) has an –OH substituent"
    # If all checks passed:
    return True, "Saturated fatty acid: contains a terminal carboxylic acid group, no C–C unsaturation, and a simple alkyl backbone"

# (Optional: main testing block)
if __name__ == "__main__":
    test_smiles = [
        "CC(=O)CCC(O)=O",   # 4-oxopentanoic acid (True)
        "CCCC(C)C(O)=O",    # 2-methylvaleric acid (True)
        "OC(C)CCCCCCCCCCCCCCCCCCC(=O)O",  # 20-hydroxyhenicosanoic acid (True)
        "CC(C)CCCCCCCCCCCCCCCCCCCCCCCCC(O)=O",  # 26-methylheptacosanoic acid (True)
        "CCCCCCCCCCCCCCCCCC(O)C([O-])=O",  # 2-hydroxyarachidate (should be rejected due to α–OH)
        "OC(=O)CCC(CCCC)CC",  # 4-Ethyloctanoic acid (should be rejected because branch is longer than methyl)
        "[O-]C(=O)CCCCCCCCC(CCCCCCCC)O",  # 10-hydroxyoctadecanoate (rejected due to extra branching)
        "OCCCCCC([O-])=O",  # 6-hydroxyhexanoate (likely rejected as not a simple fatty acid)
    ]
    for s in test_smiles:
        res, reason = is_saturated_fatty_acid(s)
        print("="*40)
        print(f"SMILES: {s}\nResult: {res}\nReason: {reason}")