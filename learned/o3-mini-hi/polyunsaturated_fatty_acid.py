"""
Classifies: CHEBI:26208 polyunsaturated fatty acid
"""
"""
Classifies: Polyunsaturated Fatty Acid
Definition: A fatty acid is defined here as a molecule having a terminal carboxylic acid group 
attached to a single, linear, non‐aromatic aliphatic chain containing more than one 
non‐aromatic C=C double bond.
Additional criteria:
  - The acid carbon in the terminal –COOH group must be attached to exactly one carbon.
  - The chain is “traced” by starting from that unique neighbor and then following the only 
    available non‐aromatic, non‐ring carbon link at each step. If at any stage a branching is encountered,
    we reject the molecule.
  - Bonds of type triple are not allowed in the chain.
  - The chain must have a minimum number of carbons (here 8) and at least two non‐aromatic C=C bonds.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyunsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule qualifies as a polyunsaturated fatty acid (PUFA).
    The criteria are:
      1. The molecule has a terminal carboxylic acid group. The carbon in the COOH must
         be bonded to exactly one carbon (thus marking the start of the chain).
      2. Starting from the chain-start carbon (the unique neighbor of the acid carbon),
         we trace a linear (non-branched), acyclic, non-aromatic chain comprising carbon atoms.
         If branching (more than one candidate) is encountered, the molecule is rejected.
      3. The chain must be at least MIN_CHAIN_LENGTH carbons long.
      4. Along the chain (i.e. between consecutive carbons in the traversed path), there must be 
         at least 2 non-aromatic double bonds (and no triple bonds).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule qualifies as a PUFA, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Identify a carboxylic acid (COOH) group using SMARTS.
    # The pattern matches a carbon (sp2) with one double-bonded oxygen and one hydroxyl.
    ca_smarts = "[CX3](=O)[O;H]"
    ca_group = Chem.MolFromSmarts(ca_smarts)
    ca_matches = mol.GetSubstructMatches(ca_group)
    if not ca_matches:
        return False, "No carboxylic acid group found; not a fatty acid"
    
    # 2. Among the acid groups, select one that is terminal.
    # That is, the acid carbon should be attached to exactly one carbon (the start of the chain).
    terminal_acid_found = False
    acid_atom = None
    chain_start = None
    for match in ca_matches:
        # match[0] is taken as the acid carbon (which is part of the COOH group)
        acid_c = mol.GetAtomWithIdx(match[0])
        # Look for neighboring carbon atoms (atomic number 6)
        carbon_neighbors = [nbr for nbr in acid_c.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) == 1:
            terminal_acid_found = True
            acid_atom = acid_c
            chain_start = carbon_neighbors[0]
            break
    if not terminal_acid_found or chain_start is None:
        return False, "Carboxylic acid group is not terminal; not a typical fatty acid"

    # 3. Trace the acyl chain starting at chain_start.
    # We “walk” along the chain by taking at each step the only possible non-aromatic,
    # non-ring carbon neighbor (excluding the atom we came from).
    chain_indices = []  # indices of carbon atoms in the chain
    current = chain_start
    previous = acid_atom  # the acid carbon is the previous atom at start
    chain_indices.append(current.GetIdx())

    while True:
        # Get neighbors that are carbons, non-aromatic, not in a ring
        candidates = []
        for nbr in current.GetNeighbors():
            # Exclude the atom we came from
            if nbr.GetIdx() == previous.GetIdx():
                continue
            if nbr.GetAtomicNum() != 6:
                continue
            if nbr.GetIsAromatic() or nbr.IsInRing():
                continue
            candidates.append(nbr)
        if len(candidates) == 0:
            # reached an endpoint
            break
        if len(candidates) > 1:
            return False, "Acyl chain is branched; expected a single linear chain for a typical fatty acid"
        # Proceed to the only candidate neighbor
        next_atom = candidates[0]
        chain_indices.append(next_atom.GetIdx())
        previous = current
        current = next_atom

    chain_length = len(chain_indices)
    MIN_CHAIN_LENGTH = 8
    if chain_length < MIN_CHAIN_LENGTH:
        return False, f"Fatty acid chain length only {chain_length} carbons; too short to be typical"
    
    # 4. Count the number of non-aromatic double bonds along the traced chain.
    double_bond_count = 0
    for i in range(len(chain_indices)-1):
        bond = mol.GetBondBetweenAtoms(chain_indices[i], chain_indices[i+1])
        if bond is None:
            continue
        # Reject if a triple bond is found in the chain.
        if bond.GetBondType() == Chem.BondType.TRIPLE:
            return False, "Chain contains a carbon–carbon triple bond which is not allowed in typical fatty acids"
        # Count non-aromatic double bonds.
        if bond.GetBondType() == Chem.BondType.DOUBLE and not bond.GetIsAromatic():
            double_bond_count += 1

    if double_bond_count <= 1:
        return False, f"Found {double_bond_count} non-aromatic C=C bond(s) in chain; need more than one to qualify as polyunsaturated"

    # Form a reason message including chain length and double bond count.
    reason = (f"Contains a terminal carboxylic acid group attached to a linear chain of {chain_length} carbons "
              f"with {double_bond_count} non-aromatic double bonds; qualifies as a polyunsaturated fatty acid")
    return True, reason


# Example usage:
if __name__ == "__main__":
    # Test with one example, resolvin D6
    test_smiles = "C(C(O)=O)C[C@@H](/C=C/C=C\\C/C=C\\C/C=C\\C=C\\[C@H](C/C=C\\CC)O)O"
    result, explanation = is_polyunsaturated_fatty_acid(test_smiles)
    print(result, explanation)