"""
Classifies: CHEBI:26607 saturated fatty acid
"""
"""
Classifies: Saturated Fatty Acid
Definition: A simple saturated fatty acid is defined (for this heuristic) as a molecule that:
  • Contains only allowed elements (C, H, O, D).
  • Is acyclic.
  • Contains exactly one carboxylic acid group (carboxyl moiety, protonated or deprotonated) whose carbon is terminal.
  • Has no carbon–carbon multiple bonds (except the carbonyl C=O in the acid group).
  • Has a simple alkyl (acyl) backbone (aside from the acid group, any branch must be a single methyl branch only;
      any –OH substituent is allowed only on the terminal (ω)–carbon).
Known adverse biological effects, when ingested to excess, are associated with these fatty acids.
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
    # Step 1. Try to parse the molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 2. Check allowed elements.
    allowed = {"C", "H", "O", "D"}
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in allowed:
            return False, f"Disallowed element found ({atom.GetSymbol()}); not a simple fatty acid"
    
    # Step 3. Reject if the molecule contains any rings.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Contains cyclic system; not a simple saturated fatty acid"
    
    # Step 4. Check that no carbon–carbon multiple bonds (double/triple) exist except in the acid group.
    # (We will later allow the C=O in the carboxyl.)
    for bond in mol.GetBonds():
        # For each bond, if it is double or triple...
        if bond.GetBondType() in [Chem.BondType.DOUBLE, Chem.BondType.TRIPLE]:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            # Allow the C=O bond inside the carboxyl group (will be handled later).
            if (a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 8) or \
               (a2.GetAtomicNum() == 6 and a1.GetAtomicNum() == 8):
                continue
            # Otherwise, if both atoms are carbons, then unsaturation exists.
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                return False, "Contains carbon–carbon unsaturation"
    
    # Step 5. Identify the unique carboxylic acid group.
    # SMARTS matches both O=C(O) and O=C([O-]) cases.
    acid_smarts = "[CX3](=O)[O;H,-]"
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if len(acid_matches) != 1:
        return False, f"Expected exactly one carboxylic acid group, found {len(acid_matches)}"
    # acid_matches returns a tuple of atom indices with the acid carbon first.
    acid_carbon_idx = acid_matches[0][0]
    acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
    
    # Step 6. Check that the acid carbon is terminal.
    carbon_neighbors = [nbr for nbr in acid_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1:
        return False, "Carboxylic acid group is not terminal (acid carbon linked to more than one carbon)"
    # The carbon neighboring the acid (α–carbon) defines the start of the acyl chain.
    alpha = carbon_neighbors[0]
    alpha_idx = alpha.GetIdx()
    
    # Step 7. Extract the main acyl chain.
    # Do a DFS starting at the α–carbon traversing through only single C–C bonds.
    def dfs_longest_chain(current_idx, visited):
        visited = visited | {current_idx}
        longest = [current_idx]
        current_atom = mol.GetAtomWithIdx(current_idx)
        for nbr in current_atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            nbr_idx = nbr.GetIdx()
            # Only traverse single bonds.
            bond = mol.GetBondBetweenAtoms(current_idx, nbr_idx)
            if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                continue
            if nbr_idx in visited:
                continue
            candidate = [current_idx] + dfs_longest_chain(nbr_idx, visited)
            if len(candidate) > len(longest):
                longest = candidate
        return longest

    main_chain = dfs_longest_chain(alpha_idx, set())
    if not main_chain:
        return False, "Could not extract acyl chain from the α–carbon"
    # For clarity, the acyl chain is the list of carbon indices from the α–carbon up to the terminal (ω) carbon.
    main_chain_set = set(main_chain)
    chain_length = len(main_chain)
    
    # Step 8. Check that the acid carbon is at one end of the full molecule.
    # We already required that the acid carbon bonds only along the chain at its α–carbon.
    # So now we check the substituents on the acyl chain.
    
    # For each carbon in the main chain, examine non-backbone substituents.
    # Allowed: at a backbone carbon, an attached carbon branch is allowed only if it is a methyl group
    #           (i.e. that branch atom should have no other heavy atom neighbors beyond the backbone attachment).
    #          an attached oxygen via a single bond is allowed only on the terminal (ω) carbon.
    for pos, atom_idx in enumerate(main_chain):
        atom = mol.GetAtomWithIdx(atom_idx)
        # Identify the two backbone neighbors in the chain (if any).
        backbone_neighbors = set()
        if pos > 0:
            backbone_neighbors.add(main_chain[pos - 1])
        if pos < chain_length - 1:
            backbone_neighbors.add(main_chain[pos + 1])
        # Also, if this is the first carbon in the chain (the α–carbon), ignore the acid carbon.
        if pos == 0:
            backbone_neighbors.add(acid_carbon_idx)
            
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in backbone_neighbors or nbr_idx in main_chain_set:
                continue
            # For a branch substitution:
            if nbr.GetAtomicNum() == 6:
                # For a methyl branch, the branch atom (nbr) must be attached only to the backbone.
                # Its degree among heavy atoms (ignoring implicit hydrogens) should be 1.
                if nbr.GetDegree() != 1:
                    return False, "Found a branch substituent longer than a methyl group on the acyl chain"
            elif nbr.GetAtomicNum() == 8:
                # Allow a single-bonded oxygen only on the terminal (ω–carbon) of the acyl chain.
                bond = mol.GetBondBetweenAtoms(atom_idx, nbr_idx)
                if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                    continue  # likely part of the acid group carbonyl; ignore
                if pos != chain_length - 1:
                    return False, "Backbone carbon (non–terminal) has a single-bonded oxygen substituent"
            else:
                return False, f"Unexpected substituent atom {nbr.GetSymbol()} on the acyl chain"

    # Step 9. (Optional additional check) Make sure the α–carbon does not have an extra –OH.
    for nbr in alpha.GetNeighbors():
        if nbr.GetAtomicNum() == 8:
            bond = mol.GetBondBetweenAtoms(alpha_idx, nbr.GetIdx())
            if bond and bond.GetBondType() == Chem.BondType.SINGLE:
                return False, "α–carbon carries a single-bonded oxygen substituent; not a typical simple saturated fatty acid"
    
    # Passed all filters.
    return True, "Saturated fatty acid: contains a terminal carboxylic acid group, no C–C unsaturation, and a simple alkyl backbone"

    
# Uncomment the block below to run some tests.
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
        "CC(C)CCCCCCCCCCCCCCCCCCCCCCCCC(O)=O",  # 26-methylheptacosanoic acid
        "OC(C)CCCCCCCCCCCCCCCCCCC(=O)O",  # 20-hydroxyhenicosanoic acid (acceptable if –OH is terminal)
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
        "CCC(O)=O",  # propionic acid
        "CCCCCCCCCCCCCCCCCCCCCCC(O)=O",  # hexacosanoic acid
        "CCCCCCCCCCCCCCCCCC(O)=O",  # hexadecanoic acid
        "OCCCCCCCCCCCCCCC(O)=O",  # 15-hydroxypentadecanoic acid
        "CCCCCCCC(O)=O",  # octanoic acid
        "CCC(C)CCCCCCCCCCCCCCCCCCCCCCC(O)=O",  # 24-methylhexacosanoic acid
        "OCCCCCCCCCCCCCCCCCCCC(O)=O",  # 20-hydroxyicosanoic acid (if the –OH is truly terminal)
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
        "CC(C)CCCCCCCCCCCCCC(O)=O",  # isotridecanoic acid
        # Examples expected to be rejected:
        "OC(=O)CCCCC(CCCCCCCCCC(CC)C)C",  # extra branch beyond methyl is present
        "OCCCCCC([O-])=O",  # OH on non–terminal carbon (not at ω–position)
        "CC(C)C(=O)CC([O-])=O",  # unsaturation outside the acid group (keto group)
    ]
    for s in test_smiles:
        res, reason = is_saturated_fatty_acid(s)
        print("="*50)
        print(f"SMILES: {s}\nResult: {res}\nReason: {reason}")