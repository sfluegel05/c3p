"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
"""
Classifies: Short-chain fatty acid.
A short-chain fatty acid is defined as an acyclic, aliphatic monocarboxylic acid that:
  – Contains exactly one carboxylic acid group (SMARTS: "[CX3](=O)[OX2H]"),
  – The acid (carboxyl) carbon is terminal (attached to exactly one other carbon),
  – All oxygen atoms present are only those in the acid group (exactly 2 oxygens),
  – And the “alkyl” chain attached to the carboxyl group (i.e. the chain starting at the unique neighbor
    of the acid carbon) must have fewer than 5 carbons connected by single bonds (so that an unbranched acid
    with 5 carbons on the alkyl chain – totaling 6) is rejected.
  – No extra non‐hydrocarbon substituents (atoms other than C or H) appear outside the carboxyl group.
  
Note: In unsaturated acids the DFS only follows single bonds so that double bonds do not “lengthen” the alkyl chain.
Some hydroxy acids or ketoacids have extra oxygens and are disqualified.
"""

from rdkit import Chem

def is_short_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule qualifies as a short-chain fatty acid.
    
    A qualifying short-chain fatty acid must:
      (a) be acyclic,
      (b) contain exactly one carboxyl group (SMARTS: "[CX3](=O)[OX2H]"),
      (c) have a terminal acid carbon (attached to exactly one other carbon),
      (d) have exactly 2 oxygens in the whole molecule,
      (e) have an alkyl chain (the saturated chain from the acid carbon’s unique carbon neighbor,
          following only single bonds) that is shorter than 5 carbons,
      (f) have no atoms outside the carboxyl functional group other than carbon or hydrogen.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as a short-chain fatty acid, False otherwise.
        str: A message explaining the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # (a) Check that the molecule is acyclic
    if mol.GetRingInfo().NumRings() != 0:
        return False, "Molecule contains rings; expected an acyclic (aliphatic) acid."
    
    # (b) Find exactly one carboxylic acid group using SMARTS "[CX3](=O)[OX2H]".
    ca_smarts = "[CX3](=O)[OX2H]"
    ca_query = Chem.MolFromSmarts(ca_smarts)
    ca_matches = mol.GetSubstructMatches(ca_query)
    if not ca_matches:
        return False, "No carboxylic acid group found."
    if len(ca_matches) != 1:
        return False, f"Expected exactly one carboxyl group, found {len(ca_matches)}."
    
    # Retrieve the matching atom indices; assume the carboxyl carbon is the first atom.
    ca_indices = ca_matches[0]
    acid_carbon_idx = ca_indices[0]
    acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
    
    # (c) Check that the acid carbon is terminal (attached to exactly one other carbon).
    carbon_neighbors = [nbr for nbr in acid_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1:
        return False, (f"The carboxyl carbon is not terminal; "
                       f"expected exactly one carbon neighbor, found {len(carbon_neighbors)}.")
    
    # (d) Count all oxygens: We require exactly 2 oxygens (from the carboxyl group).
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count != 2:
        return False, f"Expected exactly 2 oxygens (from the carboxyl group), found {oxygen_count}."
    
    # (e) Measure the alkyl chain length.
    # Only follow single bonds between carbons.
    # Start at the unique carbon neighbor attached to the acid carbon.
    start_atom_idx = carbon_neighbors[0].GetIdx()
    
    # Build a dictionary mapping atom indices (only carbons) to a list of connected carbon indices
    # but we will filter only those bonds that are SINGLE.
    carbon_only_neighbors = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        idx = atom.GetIdx()
        nbs = []
        for bond in atom.GetBonds():
            # Only consider the neighbor if the bond is a single bond.
            if bond.GetBondType() != Chem.BondType.SINGLE:
                continue
            # Get the neighbor index.
            nbr = bond.GetOtherAtom(atom)
            if nbr.GetAtomicNum() == 6:
                nbs.append(nbr.GetIdx())
        carbon_only_neighbors[idx] = nbs

    # Define DFS that only traverses carbon-only, single-bond connections.
    def dfs_max_length(current_idx, visited):
        max_len = 1  # count current atom
        for nbr in carbon_only_neighbors.get(current_idx, []):
            if nbr not in visited:
                length = 1 + dfs_max_length(nbr, visited | {nbr})
                if length > max_len:
                    max_len = length
        return max_len

    alkyl_chain_length = dfs_max_length(start_atom_idx, {start_atom_idx})
    # The chain length here excludes the acid carbon.
    if alkyl_chain_length >= 5:
        return False, (f"Molecule has an alkyl chain length of {alkyl_chain_length} (excluding the acid carbon) "
                       "which is too long for a short-chain fatty acid (unbranched analog with 6 or more carbons).")
    
    # (f) Ensure that outside the carboxyl group, no atoms other than C or H appear.
    acid_group_idxs = set(ca_indices)
    for atom in mol.GetAtoms():
        if atom.GetIdx() in acid_group_idxs:
            continue
        if atom.GetAtomicNum() not in (1, 6):
            return False, (f"Found a non-hydrocarbon substituent: atom {atom.GetSymbol()} "
                           "outside the carboxyl group is not allowed.")
    
    return True, ("Molecule contains exactly one carboxyl group with a terminal acid carbon, is acyclic, "
                  "has only the 2 oxygen atoms of the acid group, an alkyl chain (traversing only single bonds) "
                  "shorter than 5, and no extraneous heteroatoms.")

# Example usage and tests (these are taken from the provided outcome lists).
if __name__ == '__main__':
    test_cases = [
        ("CC(C)(C)C(O)=O", "pivalic acid (should be accepted)"),
        ("CCC(O)=O", "propionic acid (should be accepted)"),
        ("C=CCC(C)C(=O)O", "2-Methyl-4-pentenoic acid (should be rejected)"),
        ("OC(=O)C(CC)=CC", "2-ethyl-2-butenoic acid (should be rejected)"),
        ("OC(=O)/C=C(C)C", "3-methyl-2-enoic acid candidate (should be rejected)"),
        ("CCCC(C)C(O)=O", "2-methylvaleric acid (should be accepted)"),
        ("CCCCCC(O)=O", "n-hexanoic acid (should be rejected as unbranched)"),
        ("C(/C=C/CCC)(O)=O", "(2E)-hexenoic acid (should be accepted if unsaturation limits chain length)"),
        ("OC(C[C@H](CC)O)=O", "(R)-3-hydroxypentanoic acid (should be rejected – extra oxygen)"),
        ("OCCC(O)=O", "3-hydroxypropionic acid (should be rejected – extra oxygen)"),
        ("CCC(=O)C(O)=O", "2-oxobutanoic acid (should be rejected – extra oxygen)"),
        ("CC(C)[C@@H](C)C(O)=O", "(R)-2,3-dimethylbutyric acid (should be accepted)"),
    ]
    
    for smi, description in test_cases:
        valid, reason = is_short_chain_fatty_acid(smi)
        print(f"SMILES: {smi}\n  {description}\n  -> {valid}; Reason: {reason}\n")