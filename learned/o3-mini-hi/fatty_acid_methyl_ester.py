"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
#!/usr/bin/env python3
"""
Classifies: Fatty acid methyl ester

A fatty acid methyl ester is defined as the carboxylic ester obtained by the formal
condensation of a fatty acid with methanol. In our improved algorithm we require:
  (i) the molecule contains a methyl ester group as defined by the SMARTS "[CX3](=O)O[CH3]";
 (ii) the acyl portion – the carbon chain attached to the carbonyl (excluding the O–CH3) –
      must form a contiguous linear carbon chain of at least 7 carbons;
(iii) the chain must be “simple” meaning that any branch on the chain (aside from the connection
      to the carbonyl on the first carbon) is trivial (a single methyl, i.e. a carbon that has no
      other heavy neighbor than the chain attachment);
 (iv) the overall molecule should not be overly decorated; for a minimal fatty acid methyl ester the
      expected heavy atom count is (chain length + 4) and we allow only a small deviation.
If these conditions are met, the function returns True plus detailed reasoning.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid_methyl_ester(smiles: str):
    """
    Determines whether a molecular structure (given as SMILES) is a simple fatty acid methyl ester.
    
    The algorithm is:
     1. Identify a methyl ester group using the SMARTS "[CX3](=O)O[CH3]".
     2. For each match (the match returns indices for carbonyl carbon, carbonyl oxygen, ester oxygen, methyl carbon),
        determine the acyl carbon attached to the carbonyl (excluding the carbonyl oxygen).
     3. From that acyl carbon, perform a depth‐first search (DFS) limited to carbon atoms to find the longest
        contiguous chain. (Unsaturation is allowed since fatty acids often contain double bonds.)
     4. Check that the chain is long enough (at least 7 carbons) and “simple”. “Simple” means that at every
        atom in the chain (except that allowed connection from the first acyl carbon to the carbonyl),
        any substituent is either on the chain already or is a trivial methyl (a carbon that has no other heavy
        neighbors than the chain attachment).
     5. Finally, compare the total heavy atoms in the molecule (molecule complexity) to the minimal expected
        heavy atoms (chain_length + 4, for the carbonyl carbon, carbonyl oxygen, ester oxygen and methyl carbon)
        using a strict tolerance.
    
    Args:
      smiles (str): SMILES string for the molecule.
    
    Returns:
       (bool, str): Boolean classification and explanation.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    ester_smarts = "[CX3](=O)O[CH3]"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    if ester_pattern is None:
        return False, "Error defining ester SMARTS"

    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No methyl ester group found in the molecule"

    # Helper: recursively find the longest contiguous carbon chain
    def dfs_longest_chain(atom_idx, prev_idx, visited):
        visited.add(atom_idx)
        current_chain = [atom_idx]
        best_extension = []
        atom = mol.GetAtomWithIdx(atom_idx)
        # Only consider neighbors that are carbons and not previously visited.
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx == prev_idx:
                continue
            if nbr.GetAtomicNum() != 6:
                continue
            if nbr_idx in visited:
                continue
            chain_candidate = dfs_longest_chain(nbr_idx, atom_idx, visited.copy())
            if len(chain_candidate) > len(best_extension):
                best_extension = chain_candidate
        return current_chain + best_extension

    # Check each ester match candidate
    for match in ester_matches:
        # match order: carbonyl carbon, carbonyl O, ester O, methyl C
        carbonyl_idx, carbonylO_idx, esterO_idx, methyl_idx = match
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        # Identify the acyl neighbor: look among neighbors of the carbonyl except the carbonylO.
        acyl_neighbor = None
        for nbr in carbonyl_atom.GetNeighbors():
            if nbr.GetIdx() == carbonylO_idx:
                continue
            if nbr.GetAtomicNum() == 6:  # allow carbon (could be sp2/sp3)
                acyl_neighbor = nbr
                break
        if acyl_neighbor is None:
            continue  # try next match

        acyl_idx = acyl_neighbor.GetIdx()
        # Get the longest chain starting at the acyl neighbor. Exclude the carbonyl atom.
        chain = dfs_longest_chain(acyl_idx, carbonyl_idx, set())
        chain_length = len(chain)
        if chain_length < 7:
            # Too short to be a fatty acid chain.
            continue

        # Now check for extra decoration along the chain.
        # For each atom in the chain, look at its neighbors not in the chain.
        # For the first acyl atom, allow the carbonyl connection.
        decorated = False
        for i, a_idx in enumerate(chain):
            atom = mol.GetAtomWithIdx(a_idx)
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx in chain:
                    continue
                # For the first chain atom, allow the connection to the carbonyl.
                if i == 0 and nbr_idx == carbonyl_idx:
                    continue
                # A trivial methyl branch must be a carbon with only the bond to the chain.
                if nbr.GetAtomicNum() == 6:
                    # Count heavy neighbors (neighbors that are not hydrogen) of the branch atom excluding current atom.
                    heavy_neighbors = [a for a in nbr.GetNeighbors() if a.GetAtomicNum() > 1 and a.GetIdx() != a_idx]
                    if len(heavy_neighbors) > 0:
                        decorated = True
                        break
                else:
                    # Non-carbon substituents are not allowed.
                    decorated = True
                    break
            if decorated:
                break

        if decorated:
            # This ester match shows extra decoration beyond a simple fatty acid chain.
            continue

        # Compute extra properties to compare against a minimal fatty acid methyl ester.
        overall_rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
        mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
        actual_heavy_atoms = mol.GetNumHeavyAtoms()
        expected_heavy_atoms = chain_length + 4  # acyl chain + (carbonyl C, carbonyl O, ester O, methyl C)
        # We use a stricter tolerance to reject overly decorated molecules.
        tolerance = 2

        if abs(actual_heavy_atoms - expected_heavy_atoms) > tolerance:
            # The molecule has extra heavy atoms compared to a minimal fatty acid methyl ester.
            continue

        # Build a reason message – also reporting the chain length and rotatable bonds along the found chain.
        rotatable_in_chain = 0
        for j in range(len(chain)-1):
            bond = mol.GetBondBetweenAtoms(chain[j], chain[j+1])
            if bond and bond.GetBondTypeAsDouble() == 1.0 and not bond.IsInRing():
                rotatable_in_chain += 1

        reason = (f"Found a methyl ester group with an acyl chain of {chain_length} carbon(s); "
                  f"rotatable bonds in chain: {rotatable_in_chain}, overall rotatable bonds: {overall_rot_bonds}, "
                  f"molecule heavy atoms: {actual_heavy_atoms} (expected ~{expected_heavy_atoms} ± {tolerance}), "
                  f"molecular weight: {mol_wt:.1f} Da")
        return True, reason

    return False, "Molecule contains a methyl ester group but the acyl portion does not meet fatty acid criteria"

# Example usage:
if __name__ == "__main__":
    # A few test examples (SMILES strings taken from the provided outcomes)
    test_smiles = [
        "CCCCCC1C(C/C(=C(\\C\\C(=C(\\[H])/CCCCCCC(=O)OC)\\[H])/[H])/[H])O1",  # expected to be rejected (complex acyl portion)
        "O[C@H](CCCCCCCC(OC)=O)/C=C/C=C/CCCCC",  # expected to be rejected (has a hydroxy decoration)
        "O=C(OC)CCCCCCCCCCCC",  # methyl tridecanoate (expected positive)
        "O(C(=O)CCC(CCCC)C)C",  # formyl 4-methyl-octanoate (expected to be rejected due to extra decoration)
        "O=C(OC)CCCCCCCCCCC/C=C\\CCCCCCCC",  # methyl erucate (expected positive if chain is simple)
        "CCCCCCC(=O)OC",  # methyl heptanoate (chain too short; expected negative)
    ]
    for sm in test_smiles:
        valid, msg = is_fatty_acid_methyl_ester(sm)
        print(f"SMILES: {sm}\nResult: {valid}\nReason: {msg}\n")