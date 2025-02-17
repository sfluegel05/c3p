"""
Classifies: CHEBI:26607 saturated fatty acid
"""
"""
Classifies: Saturated Fatty Acid
Definition: A saturated fatty acid is defined for this heuristic as a molecule that:
  • Contains only allowed elements (C, H, O, D).
  • Is acyclic.
  • Contains exactly one carboxylic acid group (protonated or deprotonated) whose carbon is terminal.
  • Contains no carbon–carbon multiple bonds.
  • Has a simple alkyl (acyl) backbone – aside from the acid group,
    any carbon branch must be limited to a terminal methyl and any –OH substituent is allowed only at the terminal (ω)–carbon.
    
Known adverse biological effects, when ingested to excess, are associated with these types of fatty acids.
"""

from rdkit import Chem

def is_saturated_fatty_acid(smiles: str):
    """
    Determines whether a molecule qualifies as a simple saturated fatty acid.
    
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      bool: True if classified as a saturated fatty acid, else False.
      str: Explanation of the classification decision.
    """
    # Step 1: Parse SMILES and check allowed elements.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    # Allow only C, H, O, and D (deuterium) atoms.
    allowed = {"C", "H", "O", "D"}
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in allowed:
            return False, f"Disallowed element found ({atom.GetSymbol()}); not a simple fatty acid"
    
    # Step 2: Reject molecule if cyclic.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Contains cyclic system; not a simple saturated fatty acid"
    
    # Step 3: Identify the unique carboxylic acid group.
    # Matches either the protonated or deprotonated form.
    acid_smarts = "[CX3](=O)[O;H,-]"
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if len(acid_matches) != 1:
        return False, f"Expected exactly one carboxylic acid group, found {len(acid_matches)}"
    # acid_matches[0] returns a tuple where the first atom is the carboxyl (acid) carbon.
    acid_carbon_idx = acid_matches[0][0]
    acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
    
    # Step 4: Check that the acid carbon is terminal: it must have exactly one carbon neighbor.
    carbon_neighbors = [nbr for nbr in acid_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1:
        return False, "Carboxylic acid group is not terminal (acid carbon linked to more than one carbon)"
    alpha = carbon_neighbors[0]  # α–carbon is the first carbon of acyl chain.
    alpha_idx = alpha.GetIdx()
    
    # Step 5: Reject if any carbon–carbon multiple bonds are found (allow only C=O in the acid group).
    for bond in mol.GetBonds():
        bt = bond.GetBondType()
        if bt in [Chem.BondType.DOUBLE, Chem.BondType.TRIPLE]:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            # Skip the C=O of the acid group.
            if (a1.GetIdx() == acid_carbon_idx and a2.GetAtomicNum() == 8) or \
               (a2.GetIdx() == acid_carbon_idx and a1.GetAtomicNum() == 8):
                continue
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                return False, "Contains carbon–carbon unsaturation"
    
    # Step 6: Extract the main acyl chain.
    # Use a depth-first search (DFS) starting at the α–carbon over single C–C bonds.
    def dfs_longest_chain(current_idx, visited):
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
        return False, "Could not extract the acyl chain from the α–carbon"
    # Our DFS gives one candidate for the longest chain.
    main_chain_set = set(main_chain)
    chain_length = len(main_chain)
    
    # Step 7: Verify the main chain has only allowed substituents.
    # For each backbone carbon (in order along main_chain) check:
    #   a) Any attached carbon branch (not part of the chain or the acid group) must be a methyl only.
    #   b) Any attached oxygen (single-bonded) is allowed only on the terminal (ω)–carbon.
    for pos, atom_idx in enumerate(main_chain):
        atom = mol.GetAtomWithIdx(atom_idx)
        # Determine the backbone neighbors (previous and next in the chain).
        backbone_neighbors = set()
        if pos > 0:
            backbone_neighbors.add(main_chain[pos - 1])
        if pos < chain_length - 1:
            backbone_neighbors.add(main_chain[pos + 1])
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            # Skip if neighbor is part of the backbone or is the acid carbon (allowed attachment point).
            if nbr_idx in backbone_neighbors or nbr_idx == acid_carbon_idx or nbr_idx in main_chain_set:
                continue
            # If neighbor is carbon: it must be a methyl branch. (This branch carbon should have no other C neighbors besides the backbone attachment.)
            if nbr.GetAtomicNum() == 6:
                c_neighbors = [n.GetIdx() for n in nbr.GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIdx() != atom_idx]
                if len(c_neighbors) > 0:
                    return False, "Found a branch substituent longer than a methyl group on the acyl chain"
            # If neighbor is oxygen: allow only on the terminal (ω)–carbon.
            elif nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(atom_idx, nbr_idx)
                if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                    continue  # likely a carbonyl in the acid; ignore
                if pos != chain_length - 1:
                    return False, "Backbone carbon (non–terminal) has a single-bonded oxygen substituent"
            else:
                return False, f"Unexpected substituent atom {nbr.GetSymbol()} on the acyl chain"
                
    # Optionally, one may also check that there is no extra –OH on the α–carbon.
    for nbr in alpha.GetNeighbors():
        if nbr.GetAtomicNum() == 8:
            bond = mol.GetBondBetweenAtoms(alpha_idx, nbr.GetIdx())
            if bond and bond.GetBondType() == Chem.BondType.SINGLE:
                return False, "α–carbon carries a single-bonded oxygen substituent; not a typical simple saturated fatty acid"
    
    return True, "Saturated fatty acid: contains a terminal carboxylic acid group, no C–C unsaturation, and a simple alkyl backbone"


# Uncomment the block below to run a few tests.
if __name__ == "__main__":
    test_smiles = [
        # True positives:
        "CC(=O)CCC(O)=O",  # 4-oxopentanoic acid
        "CCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O",  # nonacosanoic acid
        "CCCC(C)C(O)=O",  # 2-methylvaleric acid
        "CC(C)CCCCCCCCCCCCCCCCCCC(O)=O",  # 20-methylhenicosanoic acid
        "C(CCCCCCCCCCCCCCCC)CCCCCCCCCCC(O)=O",  # octacosanoic acid
        "CCCCCCCCCCCCCCCCC(O)=O",  # heptadecanoic acid
        "CCCCCCCCC(C)CC(O)=O",  # 3-methylundecanoic acid
        "CC(C)C(O)=O",  # isobutyric acid
        "CCCC(O)=O",  # butyric acid
        "CC(C)CCCCCCCCCCCCCCCCCCCCC(O)=O",  # 26-methylheptacosanoic acid
        "CC(CO)CCCC(C)CCCC(C)CCCC(C)CC(O)=O",  # omega-hydroxyphytanic acid
        "CCCCCCCCCCCCC(O)=O",  # tridecanoic acid
        "CCCCCCCCCCCCCCCCCCCCCC(O)=O",  # icosanoic acid
        "OCCCCCCCCCCCCC(O)=O",  # 13-hydroxytridecanoic acid
        "CCC(C)CCCCCCCCCCCCCCCCC(O)=O",  # 18-methylicosanoic acid
        "CC(C)CCCCCCCCC(O)=O",  # 10-methylundecanoic acid
        "CC(C)CCCCCCCCCCCCCC(O)=O",  # isoheptadecanoic acid
        "CCCCCCCCCCCC(O)=O",  # dodecanoic acid
        "CC(C)CCCCCCCCCCCCCCCCCCCCC(O)=O",  # 22-methyltricosanoic acid
        "CCC(C)CCCCCCCCCCCCCCCCCCCCC(O)=O",  # 22-methyltetracosanoic acid
        "CCCCCCCCCCCCCCC[C@H](C)C(O)=O",  # (2S)-2-methylheptadecanoic acid
        "C(C(C(C(C(C(C(C(C(C(C(C(C(C(=O)O)([2H])[2H])([2H])[2H])([2H])[2H])([2H])[2H])([2H])[2H])([2H])[2H])([2H])[2H])([2H])[2H])([2H])[2H])([2H])[2H])([2H])([2H])[2H]",  # tetradecanoic-d27 acid
        "CCC(O)=O",  # propionic acid
        "CCCCCCCCCCCCCCCCCCCCCCC(O)=O",  # hexacosanoic acid
        "[2H]C([2H])([2H])C([2H])([2H])C([2H])([2H])C([2H])([2H])C([2H])([2H])C([2H])([2H])C([2H])([2H])C([2H])([2H])C([2H])([2H])C([2H])([2H])C([2H])([2H])C([2H])([2H])C([2H])([2H])C([2H])([2H])C(O)=O",  # palmitic acid-d31
        "CCCCCCCCCCCCCCCCCC(O)=O",  # hexadecanoic acid
        "OCCCCCCCCCCCCCCC(O)=O",  # 15-hydroxypentadecanoic acid
        "CCCCCCCC(O)=O",  # octanoic acid
        "CCC(C)CCCCCCCCCCCCCCCCCCCCCCC(O)=O",  # 24-methylhexacosanoic acid
        "OCCCCCCCCCCCCCCCCCCCC(O)=O",  # 20-hydroxyicosanoic acid
        "CC(C)CCCCCC(O)=O",  # 7-methyloctanoic acid
        "CC(C)C(C)C(O)=O",  # 2,3-dimethylbutyric acid
        "CC(C)CCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O",  # 28-methylnonacosanoic acid
        "CCC(C)CCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O",  # 28-methyltriacontanoic acid
        "CCCCCCCCCCCCCCCCCCCCCCCC(O)=O",  # tetracosanoic acid
        "CCCCCCC(C)CCCCCCCCCCC(O)=O",  # 12-methyloctadecanoic acid
        "CCCCCCCCCCCCC(C)C(O)=O",  # 2-methylhexadecanoic acid
        "CCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O",  # heptacosanoic acid
        "CC[C@@H](C)C(O)=O",  # (R)-2-methylbutyric acid
        "CCCCCCCCCCCCCCCCCCCCCC(O)=O",  # tricosanoic acid
        "CCC(C)C(O)=O",  # 2-methylbutyric acid
        "CC(C)(C)CC(O)=O",  # 3,3-dimethylbutyric acid
        "CCCCCC(C)CCCCCC(C)CCCCCCCCCC(C)CC(O)=O",  # 3,13,19-trimethyltricosanoic acid
        "CCC(C)CCCCCCCCCCCCC(O)=O",  # 14-methylhexadecanoic acid
        "C[C@@H](O)CCCCCCCCCCCCCCCCCCC(O)=O",  # (20R)-20-hydroxyhenicosanoic acid (may be rejected if extra OH is internal)
        "CC(C)CCCCCCCCCC(O)=O",  # isotridecanoic acid
        # Examples that should be rejected:
        "OC(=O)CCCCC(CCCCCCCCCC(CC)C)C",  # 6,16-dimethyl-octadecanoic acid (extra branch)
        "OCCCCCC([O-])=O",  # 6-hydroxyhexanoate (OH on non–terminal carbon)
        "CC(C)C(=O)CC([O-])=O",  # 4-methyl-3-oxopentanoate (contains C=O outside acid group)
        "OC(=O)CCCC(CCCCCCC)C",  # 5-methyl lauric acid (branch longer than methyl)
    ]
    for s in test_smiles:
        res, reason = is_saturated_fatty_acid(s)
        print("="*50)
        print(f"SMILES: {s}\nResult: {res}\nReason: {reason}")