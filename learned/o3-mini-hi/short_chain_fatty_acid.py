"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
"""
Classifies: Short–chain fatty acid.

A short-chain fatty acid is defined as an acyclic, aliphatic monocarboxylic acid that:
  – Contains exactly one carboxylic acid group (SMARTS: "[CX3](=O)[OX2H]"),
  – The acid (carboxyl) carbon is terminal (attached to exactly one other carbon),
  – All oxygen atoms present are only those in the acid group (exactly 2 oxygens),
  – Has a total of 3 to 6 carbon atoms,
  – And the “alkyl” chain attached to the carboxyl group (i.e. the longest continuous carbon path starting 
    from the unique carbon neighbor of the acid carbon) is less than 5 carbons long (so that adding the acid 
    carbon gives a total chain length less than 6). This disqualifies an unbranched hexanoic acid while allowing 
    a branched acid of the same total size.
  – No extra non‐hydrocarbon substituents (atoms other than C or H) appear outside the carboxyl group.
"""

from rdkit import Chem

def is_short_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule qualifies as a short-chain fatty acid.
    
    A qualifying short-chain fatty acid must:
      (a) be acyclic,
      (b) contain exactly one carboxylic acid group (SMARTS: "[CX3](=O)[OX2H]"),
      (c) have the acid carbon terminal (attached to exactly one carbon),
      (d) show that all oxygens present (count exactly 2) reside in the acid group,
      (e) contain 3–6 carbons overall,
      (f) have an alkyl chain (i.e. the longest continuous chain from the acid carbon’s unique neighbor)
          that is shorter than 5; thus an unbranched acid (with 5 carbons attached to the acid carbon, so a 
          total of 6) will be rejected,
      (g) have no non‐hydrocarbon substituents outside the carboxyl group.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule qualifies as a short-chain fatty acid, False otherwise.
        str: A message explaining the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # (a) Check that the molecule is acyclic.
    if mol.GetRingInfo().NumRings() != 0:
        return False, "Molecule contains rings; expected an acyclic (aliphatic) acid."
    
    # (b) Find exactly one carboxylic acid group using the SMARTS "[CX3](=O)[OX2H]".
    ca_smarts = "[CX3](=O)[OX2H]"
    ca_query = Chem.MolFromSmarts(ca_smarts)
    ca_matches = mol.GetSubstructMatches(ca_query)
    if not ca_matches:
        return False, "No carboxylic acid functional group found."
    if len(ca_matches) != 1:
        return False, f"Expected exactly one carboxyl group, found {len(ca_matches)}."
    
    # (c) The acid carbon is the first atom in the matching tuple.
    ca_idx = ca_matches[0][0]
    ca_atom = mol.GetAtomWithIdx(ca_idx)
    # Check that the acid carbon has exactly one carbon neighbor.
    c_neighbors = [nbr for nbr in ca_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(c_neighbors) != 1:
        return False, (f"The carboxyl carbon is not terminal; "
                       f"expected exactly one carbon neighbor, found {len(c_neighbors)}.")
    
    # (d) Count all oxygen atoms. Only 2, from the acid group, are allowed.
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count != 2:
        return False, f"Expected exactly 2 oxygens (from the carboxyl group), found {oxygen_count}."
    
    # (e) Count all carbon atoms. Must have 3–6 carbons total.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 3:
        return False, f"Total carbon count is {total_carbons}, which is too low for a fatty acid."
    if total_carbons > 6:
        return False, f"Total carbon count is {total_carbons}, which exceeds allowed short-chain length (max 6)."
    
    # (f) Assess the length of the alkyl chain attached to the carboxyl group.
    # Build a dictionary of carbon-only neighbors for every carbon atom.
    carbon_idxs = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    carbon_neighbors = {}
    for idx in carbon_idxs:
        atom = mol.GetAtomWithIdx(idx)
        # Only consider neighbors that are carbons.
        carbon_neighbors[idx] = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
    
    # Get the unique neighbor of the acid carbon.
    alkyl_start = c_neighbors[0].GetIdx()
    
    # We now compute the longest simple path (i.e. without revisiting atoms)
    # in the carbon-only graph starting from alkyl_start. Here the path length corresponds
    # to the number of carbons in the alkyl chain (excluding the acid carbon).
    def dfs_longest(current_idx, visited):
        max_length = 1  # current atom counts as length 1.
        for nbr in carbon_neighbors.get(current_idx, []):
            if nbr not in visited:
                length = 1 + dfs_longest(nbr, visited | {nbr})
                if length > max_length:
                    max_length = length
        return max_length

    alkyl_chain_length = dfs_longest(alkyl_start, {alkyl_start})
    # For an unbranched acid, the alkyl chain length equals (total carbons - 1): if it is 5 then the molecule
    # is a straight chain n-hexanoic acid which is not accepted. We require the alkyl chain length to be strictly less than 5.
    if alkyl_chain_length >= 5:
        return False, (f"Molecule has an alkyl chain length of {alkyl_chain_length} (excluding the acid carbon), "
                       "indicating an unbranched acid (n‐hexanoic acid or longer), which is not regarded as a short-chain fatty acid.")
    
    # (g) Verify that outside the carboxyl group (atoms in the matched SMARTS) only carbons and hydrogens are present.
    acid_group_idxs = set(ca_matches[0])
    for atom in mol.GetAtoms():
        if atom.GetIdx() in acid_group_idxs:
            continue  # skip atoms in the carboxyl group
        if atom.GetAtomicNum() not in (1, 6):
            return False, (f"Found a non-hydrocarbon substituent: atom {atom.GetSymbol()} "
                           "outside the carboxyl group is not allowed.")
    
    return True, ("Molecule contains one carboxyl group with a terminal acid carbon, is acyclic, "
                  "has 3–6 carbons total, an alkyl chain (from the acid carbon) shorter than 5, "
                  "and no extraneous heteroatoms.")

# Example usage:
if __name__ == '__main__':
    # These test cases are taken from the provided outcomes.
    test_cases = [
        ("CC(C)(C)C(O)=O", "pivalic acid (should be accepted)"),
        ("CCC(O)=O", "propionic acid (should be accepted)"),
        ("C=CCC(C)C(=O)O", "2-Methyl-4-pentenoic acid (should be rejected)"),
        ("OC(=O)C(CC)=CC", "2-ethyl-2-butenoic acid (should be rejected)"),
        ("OC(=O)/C=C(C)C", "3-methyl-2-enoic acid candidate (should be rejected)"),
        ("CCCC(C)C(O)=O", "2-methylvaleric acid (should be accepted)"),
        ("CCCCCC(O)=O", "n-hexanoic acid (should be rejected as unbranched)"),
        ("C(/C=C/CCC)(O)=O", "(2E)-hexenoic acid (should be accepted)"),
        ("OC(C[C@H](CC)O)=O", "(R)-3-hydroxypentanoic acid (should be rejected – extra oxygen)"),
        ("OCCC(O)=O", "3-hydroxypropionic acid (should be rejected – extra oxygen)"),
        ("CCC(=O)C(O)=O", "2-oxobutanoic acid (should be rejected – extra oxygen)"),
        ("CC(C)[C@@H](C)C(O)=O", "(R)-2,3-dimethylbutyric acid (should be accepted)"),
    ]
    
    for smi, name in test_cases:
        valid, reason = is_short_chain_fatty_acid(smi)
        print(f"SMILES: {smi}  NAME: {name}\n   -> {valid}; Reason: {reason}\n")