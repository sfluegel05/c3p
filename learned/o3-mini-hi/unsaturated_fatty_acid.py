"""
Classifies: CHEBI:27208 unsaturated fatty acid
"""
#!/usr/bin/env python
"""
Classifies: Unsaturated Fatty Acid

Definition (heuristic):
  An unsaturated fatty acid must:
    1. Be a valid molecule.
    2. Not contain any carbon–carbon triple bonds.
    3. Contain exactly one free (protonated) carboxylic acid group (matched using "C(=O)[OH]").
    4. Have its acid carbon terminal (i.e. bonded to exactly one carbon).
    5. Possess a longest acyclic (non‐ring) carbon chain starting from the only neighbor
       of the acid carbon that is at least MIN_CHAIN carbons long.
    6. Include at least one double bond along that chain.
    
Note: This heuristic method will return False for molecules with extra acid groups or where the
      chain is embedded in cycles. This additional filtering was added in response to observed false positives.
"""

from rdkit import Chem

def is_unsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule qualifies as an unsaturated fatty acid.

    The criteria are:
      1. The SMILES can be parsed.
      2. The molecule must not contain any carbon–carbon triple bonds.
      3. The molecule must contain exactly one free (protonated) carboxylic acid group,
         defined by the SMARTS "C(=O)[OH]".
      4. The acid carbon (of that group) must be terminal (connected to exactly one carbon).
      5. Starting from that single carbon neighbor, the longest acyclic (non-ring) carbon-only path
         must be at least a preset minimum length.
      6. At least one bond along that path must be a C=C double bond.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as an unsaturated fatty acid, False otherwise.
        str: Explanation of the result.
    """
    # 1. Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 2. Reject any molecule containing carbon–carbon triple bonds.
    triple_bond_pattern = Chem.MolFromSmarts("C#C")
    if mol.HasSubstructMatch(triple_bond_pattern):
        return False, "Contains carbon–carbon triple bonds which are not allowed"
    
    # 3. Find free (protonated) carboxylic acid groups using SMARTS "C(=O)[OH]".
    acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "Does not contain a free (protonated) carboxylic acid group"
    if len(acid_matches) != 1:
        return False, "Molecule should have exactly one free (protonated) carboxylic acid group"
    
    # Get the acid carbon (assumed to be the first atom in the match).
    acid_match = acid_matches[0]
    acid_c_idx = acid_match[0]
    acid_c = mol.GetAtomWithIdx(acid_c_idx)
    
    # Additionally, ensure the acid carbon is not in a ring (typical for fatty acids)
    if acid_c.IsInRing():
        return False, "Acid carbon is in a ring; not a linear fatty acid"
    
    # 4. The acid carbon must be terminal. Count its carbon neighbors.
    carbon_neighbors = [nbr.GetIdx() for nbr in acid_c.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1:
        return False, "The carboxyl group is not terminal (acid carbon is attached to more than one carbon)"
    
    # The only neighbor is the start of the fatty acid chain.
    chain_start_idx = carbon_neighbors[0]
    chain_start_atom = mol.GetAtomWithIdx(chain_start_idx)
    # If this atom is in a ring, it is not a typical fatty acyl chain.
    if chain_start_atom.IsInRing():
        return False, "The carbon chain attached to the acid group is part of a ring"

    # 5. Starting from chain_start_idx, perform a DFS that only follows carbons not in any ring.
    #    We record the chain length (number of carbons) and the bond orders along the chain.
    def dfs(current_idx, visited):
        visited.add(current_idx)
        current_atom = mol.GetAtomWithIdx(current_idx)
        # Best result: (chain_length, bond_orders along path, list of atom indices in path)
        best = (1, [], [current_idx])
        # Explore neighbors that are carbons and not in rings.
        for nbr in current_atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            # Ignore if the neighbor is in a ring.
            if nbr.IsInRing():
                continue
            nbr_idx = nbr.GetIdx()
            if nbr_idx in visited:
                continue
            bond = mol.GetBondBetweenAtoms(current_idx, nbr_idx)
            # Get the bond order (1.0 for single, 2.0 for double)
            bond_order = bond.GetBondTypeAsDouble()
            # Recursively search from this neighbor using a copy of visited.
            new_length, new_bond_orders, new_path = dfs(nbr_idx, visited.copy())
            candidate_length = 1 + new_length
            candidate_bond_orders = [bond_order] + new_bond_orders
            candidate_path = [current_idx] + new_path
            if candidate_length > best[0]:
                best = (candidate_length, candidate_bond_orders, candidate_path)
        return best
    
    MIN_CHAIN = 4  # Minimum required number of carbons in the chain (excluding the acid carbon)
    chain_length, bonds_along_chain, path = dfs(chain_start_idx, set())
    if chain_length < MIN_CHAIN:
        return False, f"Alkyl chain too short for a fatty acid (chain length = {chain_length})"
    
    # 6. Check for unsaturation. At least one bond must be a double bond (bond order 2.0).
    if not any(bond_order == 2.0 for bond_order in bonds_along_chain):
        return False, "Main carbon chain does not contain any C=C double bonds (unsaturation missing)"
    
    return True, (f"Contains a terminal carboxyl group, a sufficiently long acyclic carbon chain (length = {chain_length}), "
                  "and at least one C=C double bond in the main chain (unsaturation).")

# Example usage (for testing purposes only):
if __name__ == '__main__':
    test_smiles = [
        # True positives:
        "C(C(O)=O)C[C@@H](/C=C/C=C\\C/C=C\\C/C=C\\C=C\\[C@H](C/C=C\\CC)O)O",  # resolvin D6
        "CCCC\\C=C/C\\C=C/C=C/C=C/[C@@H]1O[C@H]1CCCC(O)=O",  # leukotriene A4
        "CCCCCCC/C=C/CCCCCCC/C(=O)O",  # (2E,4E,6E,8E,10E)-octadecapentaenoic acid-like
        "CC\\C(CC[C@H]1O[C@@]1(C)CC)=C/CC\\C(C)=C\\C(O)=O",  # juvenile hormone I acid
        "C(C(O)=O)C/C=C\\C[C@H](\\C=C\\C=C/C=C/C=C/[C@H]([C@@H](C/C=C\\CC)O)O)O",  # aspirin-triggered resolvin D2
        "CC(C)=CCCC(C)=CCCC(C)=CC(O)=O",  # farnesoic acid
        "OC(=O)CCC=C",  # pent-4-enoic acid
        "[H]C(C)=CC([H])=CC(O)=O",  # sorbic acid
        # False positive (should be rejected):
        "O=C1N(O)CCCCC1NC(=O)C(C(OC(=O)C(NC(=O)C2N=C(C3=C(O)C=CC=C3)OC2C)CCCCN(O)C(=O)/C=C\\CCCCCCCCCC(=O)O)CC)C",  # Carboxymycobactin-7
    ]
    
    for smi in test_smiles:
        result, reason = is_unsaturated_fatty_acid(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n{'-'*60}")