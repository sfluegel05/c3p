"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
"""
Classifies: Short‐chain fatty acid.

Definition:
  "An aliphatic monocarboxylic acid with a chain length of less than C6.
   If any non‐hydrocarbon substituent is present, the compound is not normally regarded as a short‐chain fatty acid."

Our implementation requires that:
  – the molecule is acyclic,
  – exactly one carboxylic acid group is present (identified by the SMARTS "[CX3](=O)[OX2H]"),
  – the carboxyl (acid) carbon is terminal in the carbon–only connectivity,
  – all oxygens in the molecule appear in that group (i.e. total oxygen count is exactly 2),
  – the total number of carbons is at least 3 and at most 6, and
  – if the acid is unbranched (i.e. the longest continuous carbon chain that includes the carboxyl carbon matches the total carbon count)
    then acids with 6 carbons (n–hexanoic acid) are rejected.
    
Note:
  In many common cases, short–chain fatty acids are branched or contain an alkene function. To discriminate them from
  unbranched acids we compute the longest path from the carboxyl carbon (in the carbon–only graph) and require that if the
  acid is completely linear then its total carbon count must be less than 6.
"""

from rdkit import Chem

def is_short_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule qualifies as a short-chain fatty acid.
    
    A qualifying short-chain fatty acid is defined as:
      (a) an acyclic, aliphatic molecule,
      (b) containing exactly one carboxylic acid group (SMARTS: "[CX3](=O)[OX2H]"),
      (c) where the acid carbon is terminal (i.e. attached to exactly one carbon in the C–only graph),
      (d) all heteroatoms (here oxygens) are only those in the COOH group,
      (e) the total number of carbons in the molecule is between 3 and 6 (inclusive),
      (f) but if the molecule is unbranched – that is, the longest continuous carbon chain (starting from the carboxyl carbon)
          covers every carbon in the molecule – then a 6–carbon acid (i.e. n–hexanoic acid) is not allowed.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as a short-chain fatty acid, False otherwise.
        str: A message describing the reason for the classification.
    """
    # Parse SMILES to molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # (a) Must be acyclic.
    if mol.GetRingInfo().NumRings() != 0:
        return False, "Molecule contains rings; expected an acyclic (aliphatic) acid."
    
    # (b) Identify the carboxylic acid group using SMARTS.
    ca_smarts = "[CX3](=O)[OX2H]"
    ca_group = Chem.MolFromSmarts(ca_smarts)
    ca_matches = mol.GetSubstructMatches(ca_group)
    if not ca_matches:
        return False, "No carboxylic acid functional group found."
    if len(ca_matches) != 1:
        return False, f"Expected exactly one carboxyl group, found {len(ca_matches)}."
    
    # For the match, by SMARTS the first atom is the carboxyl carbon.
    ca_idx = ca_matches[0][0]
    ca_atom = mol.GetAtomWithIdx(ca_idx)
    
    # (c) Check that the carboxyl carbon is terminal in the carbon-only connectivity.
    carbon_neighbors = [nbr for nbr in ca_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1:
        return False, (f"The carboxyl carbon is not terminal; "
                       f"expected exactly one carbon neighbor, found {len(carbon_neighbors)}.")
    
    # (d) Count total oxygen atoms in the molecule.
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count != 2:
        return False, f"Expected exactly 2 oxygens (from the carboxyl group), found {oxygen_count}."
    
    # (e) Count total carbon atoms.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 3:
        return False, f"Total carbon count is {total_carbons}, which is too low for a fatty acid."
    if total_carbons > 6:
        return False, f"Total carbon count is {total_carbons}, which exceeds allowed short-chain length (max 6)."
    
    # Build a list of indices for carbon atoms.
    carbon_idxs = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    
    # Create a dictionary mapping a carbon atom index to its carbon-neighbor indices.
    carbon_neighbors_dict = {}
    for idx in carbon_idxs:
        atom = mol.GetAtomWithIdx(idx)
        # Only consider neighbors that are carbon atoms.
        nbrs = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        carbon_neighbors_dict[idx] = nbrs

    # Helper: recursively compute the longest simple path in the carbon-only graph that starts at 'start_idx'
    # and does not revisit atoms.
    def dfs_longest(start_idx, visited):
        max_length = 1  # count the starting atom itself
        for nbr in carbon_neighbors_dict[start_idx]:
            if nbr not in visited:
                length = 1 + dfs_longest(nbr, visited | {nbr})
                if length > max_length:
                    max_length = length
        return max_length

    # Compute the longest path (number of carbons) that includes the carboxyl carbon.
    longest_from_ca = dfs_longest(ca_idx, {ca_idx})
    
    # (f) If the acid is completely unbranched then (by definition) total carbons equals longest chain.
    # In that case, we disallow a 6–carbon acid (i.e. n–hexanoic acid is not normally considered short chain).
    if longest_from_ca == total_carbons and total_carbons == 6:
        return False, ("Molecule is an unbranched (linear) acid with 6 carbons (n‐hexanoic acid), "
                       "which is not normally regarded as a short-chain fatty acid.")
    
    # (g) Optionally, check that outside the carboxyl group only carbons appear.
    ca_atom_idxs = set(ca_matches[0])
    for atom in mol.GetAtoms():
        if atom.GetIdx() in ca_atom_idxs:
            continue
        if atom.GetAtomicNum() not in (1, 6):  # allow hydrogen (usually implicit) and carbon
            return False, (f"Found a non-hydrocarbon substituent: atom {atom.GetSymbol()} "
                           "outside the carboxyl group is not allowed.")
    
    return True, ("Molecule contains one carboxyl group with a terminal acid carbon, is acyclic, "
                  "has 3–6 carbons total (with branched acids allowed, but unbranched 6-carbon acids rejected), "
                  "and no extraneous substituents.")

# Example usage (you can run basic tests if executing the module as main)
if __name__ == '__main__':
    # a few examples taken from the outcomes
    test_cases = [
        ("CC(C)(C)C(O)=O", "pivalic acid"),
        ("CCC(O)=O", "propionic acid"),
        ("C=CCC(C)C(=O)O", "2-Methyl-4-pentenoic acid (FP candidate)"),
        ("OC(=O)C(CC)=CC", "2-ethyl-2-butenoic acid (FP candidate)"),
        ("C/C=C/C(O)=O", "cis-pent-2-enoic acid"),
        ("CCCC(C)C(O)=O", "2-methylvaleric acid"),
        ("CCCCCC(O)=O", "n-hexanoic acid (should be rejected if unbranched)"),
    ]
    
    for smi, name in test_cases:
        valid, reason = is_short_chain_fatty_acid(smi)
        print(f"SMILES: {smi}  NAME: {name}\n   -> {valid}; Reason: {reason}\n")