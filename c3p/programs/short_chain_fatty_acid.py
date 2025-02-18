"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
"""
Classifies: Short–chain fatty acid.

Definition:
  "An aliphatic monocarboxylic acid with a chain length of less than C6.
   If any non‐hydrocarbon substituent is present, the compound is not normally regarded as a short‐chain fatty acid."

Our improved implementation requires that:
  – the molecule is acyclic,
  – exactly one carboxylic acid group is present (identified by the SMARTS "[CX3](=O)[OX2H]"),
  – the carboxyl (acid) carbon is terminal (attached to exactly one other carbon),
  – all oxygen atoms in the molecule are only those in the acid group (total oxygen count is exactly 2),
  – the total number of carbons is at least 3 and at most 6,
  – and the “chain length” (the longest continuous path in the carbon–only graph starting from the carboxyl carbon)
    is less than 6. (This way, an unbranched 6–carbon acid is rejected while a branched acid may be allowed.)
    
If any extra non‐hydrocarbon (i.e. non C/H) substituent is found outside the carboxyl group, the molecule is rejected.
"""

from rdkit import Chem

def is_short_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule qualifies as a short-chain fatty acid.
    
    A qualifying short-chain fatty acid is defined as:
      (a) an acyclic, aliphatic molecule,
      (b) containing exactly one carboxylic acid group (SMARTS: "[CX3](=O)[OX2H]"),
      (c) where the acid carbon is terminal (i.e. attached to exactly one carbon),
      (d) all oxygens in the molecule are in the carboxyl group (so total oxygen count must be 2),
      (e) the molecule has 3–6 carbons in total,
      (f) and the longest continuous carbon chain starting at the acid carbon is less than 6 (i.e. an unbranched hexanoic acid is disallowed).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as a short-chain fatty acid, False otherwise.
        str: A message describing the reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # (a) The molecule must be acyclic.
    if mol.GetRingInfo().NumRings() != 0:
        return False, "Molecule contains rings; expected an acyclic (aliphatic) acid."
    
    # (b) Locate exactly one carboxylic acid group using SMARTS "[CX3](=O)[OX2H]"
    ca_smarts = "[CX3](=O)[OX2H]"
    ca_group = Chem.MolFromSmarts(ca_smarts)
    ca_matches = mol.GetSubstructMatches(ca_group)
    if not ca_matches:
        return False, "No carboxylic acid functional group found."
    if len(ca_matches) != 1:
        return False, f"Expected exactly one carboxyl group, found {len(ca_matches)}."
    
    # The first atom in the match is the carboxyl carbon.
    ca_idx = ca_matches[0][0]
    ca_atom = mol.GetAtomWithIdx(ca_idx)
    
    # (c) Verify that the acid carbon is terminal in the carbon-only connectivity.
    # It must have exactly one neighbor which is a carbon.
    c_neighbors = [nbr for nbr in ca_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(c_neighbors) != 1:
        return False, (f"The carboxyl carbon is not terminal; expected exactly one carbon neighbor, "
                       f"found {len(c_neighbors)}.")
    
    # (d) Count all oxygen atoms in the molecule.
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count != 2:
        return False, f"Expected exactly 2 oxygens (from the carboxyl group), found {oxygen_count}."
    
    # (e) Count all carbon atoms.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 3:
        return False, f"Total carbon count is {total_carbons}, which is too low for a fatty acid."
    if total_carbons > 6:
        return False, f"Total carbon count is {total_carbons}, which exceeds allowed short-chain length (max 6)."
    
    # (f) Build the carbon-only graph (mapping from carbon index to carbon neighbor indices).
    carbon_idxs = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    carbon_neighbors = {}
    for idx in carbon_idxs:
        atom = mol.GetAtomWithIdx(idx)
        # Only include neighbors that are carbon; ignore oxygens and others.
        carbon_neighbors[idx] = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
    
    # Helper function: recursively compute the longest simple (acyclic) path in the carbon graph
    # starting from a given carbon using depth-first search.
    def dfs_longest(current_idx, visited):
        max_length = 1
        for nbr in carbon_neighbors[current_idx]:
            if nbr not in visited:
                length = 1 + dfs_longest(nbr, visited | {nbr})
                if length > max_length:
                    max_length = length
        return max_length

    longest_chain = dfs_longest(ca_idx, {ca_idx})
    # For a short-chain acid, even if the molecule is branched (total carbons ≤ 6),
    # the longest continuous chain starting at the acid carbon must be less than 6.
    if longest_chain >= 6:
        return False, ("Molecule is unbranched with a continuous carbon chain length "
                       f"of {longest_chain} (n‐hexanoic acid or longer), which is not regarded as a short-chain fatty acid.")
    
    # (g) Verify that outside the carboxyl group (the atoms in the matched SMARTS) only carbons and hydrogens appear.
    # (This excludes any extraneous non‐hydrocarbon substituents.)
    ca_atom_idxs = set(ca_matches[0])
    for atom in mol.GetAtoms():
        if atom.GetIdx() in ca_atom_idxs:
            continue
        # Allow only carbon (6) or hydrogen (1); note that implicit hydrogens are not present in the atom list.
        if atom.GetAtomicNum() not in (1, 6):
            return False, (f"Found a non-hydrocarbon substituent: atom {atom.GetSymbol()} "
                           "outside the carboxyl group is not allowed.")
    
    return True, ("Molecule contains one carboxyl group with a terminal acid carbon, is acyclic, "
                  "has 3–6 carbons total, and a continuous carbon chain (from the acid carbon) shorter than 6; "
                  "there are no extraneous heteroatoms.")

# Example usage when executing the module as main.
if __name__ == '__main__':
    # Sample test cases (names taken from the outcomes provided)
    test_cases = [
        ("CC(C)(C)C(O)=O", "pivalic acid"),
        ("CCC(O)=O", "propionic acid"),
        ("C=CCC(C)C(=O)O", "2-Methyl-4-pentenoic acid (false positive candidate)"),
        ("OC(=O)C(CC)=CC", "2-ethyl-2-butenoic acid (false positive candidate)"),
        ("C/C=C/C(O)=O", "cis-pent-2-enoic acid"),
        ("CCCC(C)C(O)=O", "2-methylvaleric acid"),
        ("CCCCCC(O)=O", "n-hexanoic acid (should be rejected if unbranched)"),
        ("OC(C[C@H](CC)O)=O", "(R)-3-hydroxypentanoic acid (should be rejected – extra oxygen)"),
    ]
    
    for smi, name in test_cases:
        valid, reason = is_short_chain_fatty_acid(smi)
        print(f"SMILES: {smi}  NAME: {name}\n   -> {valid}; Reason: {reason}\n")