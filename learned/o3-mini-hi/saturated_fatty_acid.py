"""
Classifies: CHEBI:26607 saturated fatty acid
"""
"""
Classifies: Saturated Fatty Acid
Definition: A saturated fatty acid (for this heuristic) must:
  1. Have a valid SMILES made only of C, H, O and D.
  2. Be acyclic.
  3. Contain exactly one terminal carboxylic acid group 
     (matched by C(=O)[O;H,-]) where the acid carbon is attached to exactly one other carbon.
  4. Have no carbon–carbon multiple bonds (C=C or C≡C).
  5. Have an “alkyl” (acyl) backbone that is “simple”, meaning that 
     aside from the acid group at one end, any branch off the main chain must be only a methyl group.
  6. Allow a single-bonded oxygen (–OH) only on the terminal (ω)–carbon.
Note: Some molecules that are “modified” (e.g. additional –OH on an internal carbon, or branches longer than methyl)
      are rejected.
      
The approach:
  • Verify allowed elements and acyclicity.
  • Locate the unique carboxylic acid group via SMARTS that accepts either protonated or deprotonated oxygen.
  • Confirm that the acid carbon is terminal (exactly one C neighbor).
  • Starting from the α–carbon, perform a DFS (depth-first search) to extract the longest carbon chain.
  • For every backbone carbon, check that any attached carbon branch (i.e. not in the main chain) is a methyl only
    and that any single-bonded oxygen substituent is found only on the last (ω)–carbon.
  • Look for any non–allowed C=C or C≡C bonds.
  
If any check fails, a message is returned stating the reason.
"""

from rdkit import Chem

def is_saturated_fatty_acid(smiles: str):
    """
    Determines whether a molecule qualifies as a simple saturated fatty acid.
    
    Args:
      smiles (str): SMILES string of the molecule
      
    Returns:
      bool: True if classified as a saturated fatty acid, else False.
      str: Explanation of the classification decision.
    """
    # Step 1: Parse SMILES, check allowed elements (only C, H, O, D)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    allowed = {"C", "H", "O", "D"}
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in allowed:
            return False, f"Disallowed element found ({atom.GetSymbol()}); not a simple fatty acid"
    
    # Step 2: Reject if molecule is cyclic.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Contains cyclic system; not a simple saturated fatty acid"
    
    # Step 3: Identify the carboxylic acid group.
    # Using SMARTS that matches either the protonated or deprotonated form.
    acid_smarts = "[CX3](=O)[O;H,-]"
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if len(acid_matches) != 1:
        return False, f"Expected exactly one carboxylic acid group, found {len(acid_matches)}"
    # acid_matches[0]: first atom is the acid carbon.
    acid_match = acid_matches[0]
    acid_carbon_idx = acid_match[0]
    acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
    # Step 4: Check that the acid carbon is terminal (has exactly one carbon neighbor)
    carbon_neighbors = [nbr for nbr in acid_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1:
        return False, "Carboxylic acid group is not terminal (acid carbon is linked to more than one carbon)"
    alpha = carbon_neighbors[0]  # α–carbon (first carbon of the acyl chain)
    alpha_idx = alpha.GetIdx()
    
    # Rule: α–carbon should not be substituted by a singly bonded oxygen (i.e. an –OH) 
    # because that would violate the “simple backbone” rule.
    for nbr in alpha.GetNeighbors():
        if nbr.GetAtomicNum() == 8:
            bond = mol.GetBondBetweenAtoms(alpha_idx, nbr.GetIdx())
            if bond and bond.GetBondType() == Chem.BondType.SINGLE:
                return False, "α–carbon carries a single-bonded oxygen substituent; not a typical saturated fatty acid"
    
    # Step 5: Reject if any C=C or C≡C bonds (except C=O in the acid group)
    for bond in mol.GetBonds():
        bt = bond.GetBondType()
        if bt in [Chem.BondType.DOUBLE, Chem.BondType.TRIPLE]:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                # Allow if one is in the acid carbon and its partner is oxygen via a double bond? (Not possible here)
                return False, "Contains carbon–carbon unsaturation"
    
    # Step 6: Extract the main acyl chain.
    # We perform a DFS from the α–carbon over C–C bonds.
    def dfs_longest_chain(current_idx, visited):
        # Returns the longest chain (list of carbon indices) starting at current_idx.
        visited = visited | {current_idx}
        longest = [current_idx]
        current_atom = mol.GetAtomWithIdx(current_idx)
        for nbr in current_atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            nbr_idx = nbr.GetIdx()
            if nbr_idx in visited:
                continue
            candidate = [current_idx] + dfs_longest_chain(nbr_idx, visited)
            if len(candidate) > len(longest):
                longest = candidate
        return longest

    main_chain = dfs_longest_chain(alpha_idx, set())
    if len(main_chain) < 1:
        return False, "Could not extract an acyl carbon chain from the α–carbon"
    
    # For debugging, sort the main chain in order (our DFS gives one candidate longest path).
    # Our main chain is assumed to run from the acid-attached carbon (α) to the ω–end.
    # (We assume that the chosen longest path is the acyl backbone.)
    # Now, create a set for easy membership checks.
    main_chain_set = set(main_chain)
    
    # Step 7: Check substitution on the main chain.
    # For each carbon in the backbone, look for substituents that are not in the chain (and also not the acid carbon).
    # Allowed: if a neighbor is a carbon then it must be a terminal methyl (i.e. it has no other carbon neighbor besides the backbone).
    # Also, a single-bonded oxygen substituent is allowed only at the terminal (ω)–carbon.
    chain_length = len(main_chain)
    for i, atom_idx in enumerate(main_chain):
        atom = mol.GetAtomWithIdx(atom_idx)
        # Get set of neighbors that are in the main chain (previous and next will normally be there)
        backbone_neighbors = set()
        if i > 0:
            backbone_neighbors.add(main_chain[i-1])
        if i < chain_length - 1:
            backbone_neighbors.add(main_chain[i+1])
        
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            # Skip if neighbor is acid carbon (only valid as the terminal acid attachment)
            if nbr_idx == acid_carbon_idx:
                continue
            # Skip if neighbor is in the validated backbone
            if nbr_idx in backbone_neighbors or nbr_idx in main_chain_set:
                continue
            # If neighbor is carbon, check that it is just a single methyl group (i.e. it is terminal in the C–C subgraph)
            if nbr.GetAtomicNum() == 6:
                # Count how many carbon neighbors nbr has (excluding the one connecting to the backbone)
                carbon_neighbors = [n.GetIdx() for n in nbr.GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIdx() != atom_idx]
                if len(carbon_neighbors) != 0:
                    return False, "Found a branch substituent longer than a methyl group on the acyl chain"
            # If neighbor is oxygen, then the bond must be a double bond (as in the acid carbonyl),
            # or—if single-bonded—the oxygen is allowed only on the terminal (ω)–carbon.
            elif nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(atom_idx, nbr_idx)
                if bond.GetBondType() != Chem.BondType.DOUBLE:
                    if i != chain_length - 1:
                        return False, "Backbone carbon (non–terminal) has a single-bonded oxygen substituent"
    # All checks passed.
    return True, "Saturated fatty acid: contains a terminal carboxylic acid group, no C–C unsaturation, and a simple alkyl backbone"

# (Optional main block for testing a few examples)
if __name__ == "__main__":
    test_smiles = [
        # True positives (examples)
        "CC(=O)CCC(O)=O",   # 4-oxopentanoic acid
        "CCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O",   # nonacosanoic acid
        "CCCC(C)C(O)=O",    # 2-methylvaleric acid
        "CC(C)CCCCCCCCCCCCCCCCCCC(O)=O",  # 20-methylhenicosanoic acid
        "C(CCCCCCCCCCCCCCCC)CCCCCCCCCCC(O)=O",  # octacosanoic acid
        "CCCCCCCCCCCCCCCCC(O)=O",  # heptadecanoic acid
        "CCCCCCCCC(C)CC(O)=O",  # 3-methylundecanoic acid
        "CC(C)C(O)=O",  # isobutyric acid
        "CCCC(O)=O",   # butyric acid
        # A couple of examples that were false positives in the earlier attempt (should fail now):
        "OC(=O)CCCCC(CCCCCCCCCC(CC)C)C",  # 6,16-dimethyl-octadecanoic acid
        "OCCCCCC([O-])=O",  # 6-hydroxyhexanoate
    ]
    for s in test_smiles:
        res, reason = is_saturated_fatty_acid(s)
        print("="*50)
        print(f"SMILES: {s}\nResult: {res}\nReason: {reason}")