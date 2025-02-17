"""
Classifies: CHEBI:39418 straight-chain saturated fatty acid
"""
#!/usr/bin/env python
"""
Classifies: Straight-chain saturated fatty acid
Definition: Any saturated fatty acid lacking a side-chain.
Here we require that the molecule has exactly one carboxyl group (–COOH) and that the carbon-only
skeleton (i.e. ignoring non-carbon atoms) forms one linear chain. In addition, extra heteroatom substituents 
attached at a non-terminal carbon (except for a carbonyl group—that is, a C=O) are not allowed.
"""

from rdkit import Chem

def is_straight_chain_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a straight-chain saturated fatty acid.
    (A) The molecule must contain exactly one carboxyl group ([CX3](=O)[OX1H0,OX2H1]).
    (B) All carbon atoms connected by single bonds should form one simple, unbranched chain.
    (C) Extra substituents on the main chain are not allowed except:
          - For the carboxyl carbon: only the two oxygens of the carboxyl group.
          - For an end (methyl) carbon: an –OH group is tolerated.
          - For an internal carbon: only a double-bonded oxygen (C=O) is tolerated.
    (D) No carbon–carbon bonds may be unsaturated.
    
    Args:
         smiles (str): a SMILES string representing the molecule.
         
    Returns:
         (bool, str): True plus an explanatory message if the molecule is classified as a 
         straight-chain saturated fatty acid; otherwise, False and a reason message.
    """
    # Parse molecule from SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # (D) Check that there are no unsaturated carbon-carbon bonds.
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
            # Only single bonds are allowed in a saturated fatty acid.
            if bond.GetBondType() != Chem.BondType.SINGLE:
                return False, "Found an unsaturated (non-single) carbon–carbon bond"
    
    # (A) Find the carboxyl group.
    # This SMARTS finds a carbon that is double bonded to an oxygen and single-bonded to an oxygen (–OH/–O–)
    carboxyl_smarts = "[CX3](=O)[OX1H0,OX2H1]"
    carboxyl_pattern = Chem.MolFromSmarts(carboxyl_smarts)
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if len(carboxyl_matches) != 1:
        return False, f"Found {len(carboxyl_matches)} carboxyl groups, expected exactly 1"
    # Use the first match. By convention, we take the first atom (the carbon) as the carboxyl carbon.
    carboxyl_match = carboxyl_matches[0]
    carboxyl_carbon_idx = carboxyl_match[0]
    # Save the allowed neighbors for the carboxyl carbon (the two oxygens in the COOH group)
    allowed_for_carboxyl = set(carboxyl_match[1:])  # typically the carbonyl O and hydroxyl O

    # (B) Build a connectivity graph of carbon atoms (atomic number 6) using only single bonds.
    carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    # Create a dictionary mapping each carbon index to its carbon neighbors.
    carbon_graph = {idx: [] for idx in carbon_indices}
    for bond in mol.GetBonds():
        # Consider only single bonds between two carbons.
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6 and bond.GetBondType() == Chem.BondType.SINGLE:
            idx1 = a1.GetIdx()
            idx2 = a2.GetIdx()
            if idx1 in carbon_graph and idx2 in carbon_graph:
                carbon_graph[idx1].append(idx2)
                carbon_graph[idx2].append(idx1)
                
    if carboxyl_carbon_idx not in carbon_graph:
        return False, "Carboxyl carbon not found in carbon connectivity graph"
    
    # Helper: DFS to get longest simple path starting from start_node.
    def dfs_longest_path(current, visited):
        longest = visited
        for nbr in carbon_graph[current]:
            if nbr not in visited:
                new_path = dfs_longest_path(nbr, visited + [nbr])
                if len(new_path) > len(longest):
                    longest = new_path
        return longest
    
    # Find the longest path starting from the carboxyl carbon.
    main_chain = dfs_longest_path(carboxyl_carbon_idx, [carboxyl_carbon_idx])
    
    # (B cont.) For a straight chain fatty acid, all carbon atoms should lie on the main chain.
    if set(main_chain) != set(carbon_graph.keys()):
        return False, "Extra carbon atoms (branch(es)) detected outside the main chain"
    
    # (C) Now check each carbon in the main chain for improper substituents.
    # First, get a mapping from atom index to the corresponding Atom object.
    idx2atom = {atom.GetIdx(): atom for atom in mol.GetAtoms()}
    # Get the list of carbon neighbors for each carbon (including non-carbon substituents).
    # We use the full molecular connectivity here.
    chain_length = len(main_chain)
    for pos, c_idx in enumerate(main_chain):
        atom = idx2atom[c_idx]
        # Determine the “role” of this chain carbon.
        is_carboxyl = (c_idx == carboxyl_carbon_idx)
        # Terminal if it is first or last in the chain.
        is_terminal = (pos == 0 or pos == chain_length - 1)
        # Examine all neighbors (from the full molecule) that are not part of the main chain.
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in main_chain:
                continue  # this bond is part of the main chain
            # Get the bond between atom and nbr.
            bond = mol.GetBondBetweenAtoms(c_idx, nbr_idx)
            # For the carboxyl carbon, only allow the oxygens in the carboxyl group.
            if is_carboxyl:
                if nbr_idx not in allowed_for_carboxyl:
                    return False, f"Carboxyl carbon has an extra substituent (atom {nbr_idx})"
                else:
                    continue
            # For an internal (nonterminal) carbon, allow one extra substituent only if it is a double-bonded oxygen (a keto group).
            # (Note: single-bonded heteroatoms on an internal carbon are considered a side‐chain.)
            if not is_terminal:
                if bond.GetBondType() == Chem.BondType.DOUBLE and nbr.GetAtomicNum() == 8:
                    continue  # allowed keto group
                else:
                    return False, f"Internal chain carbon (atom {c_idx}) has a disallowed substituent (atom {nbr_idx})"
            else:
                # For the terminal (non-carboxyl) end, typically a CH3 is expected.
                # Allow an extra single-bonded oxygen (e.g. –OH) at the terminal end (omega–hydroxy fatty acid),
                # but do not allow more than one such substituent.
                # Count the number of nonchain substituents.
                # (We already are inside a loop over such neighbors.)
                if nbr.GetAtomicNum() == 8:
                    # If the bond is double, that is a keto; if single, that is an -OH.
                    # For a terminal methyl, a single -OH is allowed.
                    if bond.GetBondType() in (Chem.BondType.SINGLE, Chem.BondType.DOUBLE):
                        # allowed (if there is only one such substituent, the loop will check them individually)
                        continue
                    else:
                        return False, f"Terminal carbon (atom {c_idx}) has an improper bond to oxygen (atom {nbr_idx})"
                else:
                    return False, f"Terminal carbon (atom {c_idx}) has an unexpected substituent (atom {nbr_idx}, atomic number {nbr.GetAtomicNum()})"

    return True, "Molecule is a straight-chain saturated fatty acid with no side-chain substituents"

# Example usage: (for testing purposes)
if __name__ == "__main__":
    test_smiles = [
        # True positives:
        "CC(=O)CCC(O)=O",                        # 4-oxopentanoic acid (correct)
        "CCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O",        # heptocasanoic acid (correct)
        "CCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O",       # nonocasanoic acid (correct)
        "CCCCCCCCCCCCCCCCCCCCCC(O)=O",            # tricosanoic acid (correct)
        # isotopically labeled (but structure identical in connectivity):
        "C(C(C(C(C(C(C(C(C(C(C(C(C(C(=O)O)([2H])[2H])([2H])[2H])([2H])[2H])([2H])[2H])([2H])[2H])([2H])[2H])([2H])[2H])([2H])[2H])([2H])[2H])([2H])[2H])",  # tetradecanoic-d27 acid
        "C(CCCCCCCCCCCCCCCC)CCCCCCCCCCC(O)=O",    # octacosanoic acid (correct)
        "CCCCCCCCCCCCCCCC(O)=O",                  # heptadecanoic acid (correct)
        "CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O",      # hexacosanoic acid (correct)
        "OC(C)CCCCCCCCCCCCCCCCCCC(O)=O",           # (20R)-20-hydroxyhenicosanoic acid (correct)
        "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O",     # triacontanoic acid (correct)
        "[2H]C([2H])([2H])C([2H])([2H])C([2H])([2H])C([2H])([2H])C([2H])([2H])C([2H])([2H])C([2H])([2H])C([2H])([2H])C(O)=O",  # palmitic acid-d31 (correct)
        "CCCCCCCCCCCCCCCC(O)=O",                  # hexadecanoic acid (correct)
        "CCCCCCC(O)=O",                          # heptanoic acid (correct)
        "OC(C)CCCCCCCCCCCCCCCCCCC(=O)O",           # 20-hydroxyhenicosanoic acid (another representation, correct)
        "OCCCCCCCCCCCCCCC(O)=O",                  # 15-hydroxypentadecanoic acid (correct)
        "CCCC(O)=O",                             # butyric acid (correct)
        "CCCCCCCCCCCCCCCCCCCC(O)=O",              # henicosanoic acid (correct)
        "CCCCCCCCCCCCCCCCCCCCCCCCC(O)=O",         # docosanoic acid (correct)
        "CCCCCCCC(O)=O",                         # octanoic acid (correct)
        "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O",    # dotriacontanoic acid (correct)
        "CCCCCCCCCCCCCC(O)=O",                    # tetradecanoic acid (correct)
        "OCCCCCCCCCCCCCCCCCCCC(O)=O",             # 20-hydroxyicosanoic acid (correct)
        "CCCCCCCCCCCCCCCCCCC(O)=O",               # nonadecanoic acid (correct)
        "CCCCCC(O)=O",                           # hexanoic acid (correct)
        # False positives (should fail):
        "CCCCCCCCCCCCCCCCCC(O)C([O-])=O",          # 2-hydroxyarachidate (should be rejected because the OH is on an internal carbon)
        "[O-]C(=O)CCCCCCCCC(CCCCCCCC)O",           # 10-hydroxyoctadecanoate (branched carbon chain)
        "CCCCCCCC[C@@H](O)[C@@H](O)CCCCCCCC([O-])=O",  # dihydroxyoctadecanoate (multiple hydroxy groups internally)
        "O[C@H](CCCCCC)C(O)=O",                    # (R)-2-Hydroxycaprylic acid (OH at C2 disqualifies)
        "O(P([O-])([O-])=O)CC([C@H](C(=O)[O-])O)=O", # (R)-2-hydroxy-3-oxo-4-(phosphonatooxy)butanoate (extra substituents)
        "CCCC(N)C(O)=O",                          # 2-aminopentanoic acid (extra amino substituent)
        "OC[C@@H](O)[C@@H](O)[C@@H](O)C(=O)C([O-])=O", # beta-D-galacturonate (sugar acid, branched)
        "C[C@H](O)CC([O-])=O",                     # (S)-3-hydroxybutyrate (OH on internal carbon)
        "OCCCCCC([O-])=O",                        # 6-hydroxyhexanoate (OH on non-terminal position)
        "NC(CCS)C(O)=O",                          # homocysteine (extra amino and thiol substituents)
        "C(CCCCCCCCCCCCCC)C[C@H](C([O-])=O)OO",    # (2R)-2-hydroperoxyoctadecanoate (branching hydroperoxy)
        "OC[C@@H](O)[C@H](O)[C@H](O)C(=O)C([O-])=O", # 2-dehydro-3-deoxy-D-gluconic acid (branched sugar acid)
        "OC(=O)C1CCC(=O)O1",                       # 5-oxotetrahydrofuran-2-carboxylic acid (ring structure)
        "OC(=O)C(N)CCCCCCCCCCCC",                  # 2-Amino-tetradecanoic acid (amino substituent on an internal carbon)
        "OP(=O)(OC[C@]1([C@H]([C@@H](CC(O1)C(=O)O)O)O)O)O",  # (2xi)-3-deoxy-7-O-phospho-beta-D-threo-hept-6-ulopyranosonic acid (extra phosphate)
        "[C@@H](N)(CN)C(=O)O",                     # 3-amino-D-alanine (simple amino acid, not a fatty acid)
        "OC(=O)C[C@H](N)CCCCCC",                   # 3R-Aminononanoic acid (amino substituent)
        "OC(=O)C(F)F",                            # difluoroacetic acid (extra halogens)
        "C[C@@H](N)C(O)=O",                        # D-alanine (amino acid)
        "OCCCCCCCCCCCCCCCCCCCCCCCCCC([O-])=O",      # 26-hydroxyhexacosanoate (OH on internal carbon)
        "OC(=O)CC(N)CCCCN",                        # 3,7-Diaminoheptanoic Acid (extra amino groups)
        "[O-]C(=O)[C@H]1CCC(=O)N1",                # 5-oxo-D-prolinate (ring-containing, extra substituents)
        "OC[C@@H](O)[C@H](O)[C@H](O)C(=O)C([O-])=O", # 2-dehydro-D-galactonate (sugar acid, branched)
        "OC(CCCCCC(O)=O)CCC",                      # 7-hydroxy-decanoic acid (branch on internal carbon)
        "OCCCCCCCCCCCCCCCCCCCCCC(O)=O",            # 22-hydroxydocosanoic acid (OH on internal carbon)
        # False negative example (should be rejected)
        "C(C[N+](C)(C)C)OP(=O)([O-])OC[C@H](OC)COC(CCCCCCCCCCCCCCC)=O", # 1-hexadecanoyl-2-methyl-sn-glycero-3-phosphocholine (no free carboxyl)
    ]
    for smi in test_smiles:
        valid, reason = is_straight_chain_saturated_fatty_acid(smi)
        print(f"SMILES: {smi}\n  -> {valid}: {reason}\n")