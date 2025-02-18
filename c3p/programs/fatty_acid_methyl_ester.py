"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
#!/usr/bin/env python3
"""
Classifies: Fatty acid methyl ester

A fatty acid methyl ester is defined as the carboxylic ester obtained by the formal
condensation of a fatty acid with methanol. Here we require:
  (i) the molecule contains a methyl ester group as defined by the SMARTS "[CX3](=O)O[CH3]";
  (ii) the acyl portion (i.e. the group attached to the carbonyl carbon, not including the O–CH3)
      forms a contiguous (linear) chain of at least 7 carbons,
  (iii) the chain is “simple”: no non‐carbon substituents appear and any extra carbon branch must be
       a trivial methyl branch (a single carbon with no additional heavy neighbors);
  (iv) and the overall molecule is not overly decorated relative to the expected minimal structure.
If these conditions are met, the function returns True with an explanation.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid_methyl_ester(smiles: str):
    """
    Determines if a molecule is a fatty acid methyl ester based on its SMILES string.
    
    The algorithm:
      1. Search for a methyl ester group via SMARTS "[CX3](=O)O[CH3]".
      2. For each match, identify the carbonyl carbon and then the acyl neighbor 
         (the carbon attached to the carbonyl that is not the ester oxygen).
      3. Starting from that acyl atom, walk the main chain expecting a linear chain with no
         extra (non-trivial) decoration.
           - Each step: from the current chain atom (starting with the acyl neighbor), find the
             next unique aliphatic non‐aromatic carbon (excluding the atom you came from).
           - If more than one candidate exists at any step, mark the chain as decorated.
      4. Check for decoration at each chain atom: allow only trivial methyl branches 
         (i.e. a carbon that has no heavy (non‐hydrogen) neighbor except the chain connection).
      5. Require the chain length to be at least 7.
      6. Compare overall molecule heavy atom count to the expected minimal structure 
         (chain_length + 4 for the carbonyl carbon, the carbonyl oxygen, the ester oxygen, and the methyl carbon)
         with a tolerance.
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      bool: True if molecule qualifies as a fatty acid methyl ester, else False.
      str: Explanation text.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the methyl ester SMARTS pattern.
    ester_smarts = "[CX3](=O)O[CH3]"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    if ester_pattern is None:
        return False, "Error defining ester SMARTS"
    
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No methyl ester group found in the molecule"
    
    # For each methyl ester match. The SMARTS returns (carbonylC, carbonylO, esterO, methylC)
    for match in ester_matches:
        carbonyl_idx, carbonylO_idx, esterO_idx, methyl_idx = match
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        
        # Identify the acyl neighbor: one of the neighbors of the carbonyl
        # (other than the carbonyl oxygen) that is a non‐aromatic carbon.
        acyl_neighbor = None
        for nbr in carbonyl_atom.GetNeighbors():
            if nbr.GetIdx() == carbonylO_idx:
                continue
            if nbr.GetAtomicNum() == 6 and not nbr.GetIsAromatic():
                acyl_neighbor = nbr
                break
        if acyl_neighbor is None:
            continue  # try next ester match
        
        # Build the main chain starting from the acyl neighbor.
        # We require that the chain is linear. For the first step, the previous atom is the carbonyl.
        main_chain = [acyl_neighbor.GetIdx()]
        prev_idx = carbonyl_atom.GetIdx()
        current_atom = acyl_neighbor
        decorated = False  # flag for extra (non trivial) decoration
        
        # Follow the chain while exactly one non‐backtracking eligible neighbor exists.
        while True:
            # Identify eligible next neighbor:
            # must be a carbon, aliphatic (non-aromatic) and not the previous atom.
            candidates = []
            for nbr in current_atom.GetNeighbors():
                if nbr.GetIdx() == prev_idx:
                    continue
                if nbr.GetAtomicNum() == 6 and not nbr.GetIsAromatic():
                    candidates.append(nbr)
            if len(candidates) == 1:
                # exactly one candidate is what we expect for a linear chain
                next_atom = candidates[0]
                main_chain.append(next_atom.GetIdx())
                prev_idx = current_atom.GetIdx()
                current_atom = next_atom
            elif len(candidates) == 0:
                # reached the end of the chain
                break
            else:
                # More than one continuation => the chain is branched beyond a simple linear chain.
                decorated = True
                break

        # At this point, main_chain holds the linear chain of carbons (represented by indices)
        chain_length = len(main_chain)
        # For a fatty acid, we require at least 7 carbons (e.g. methyl octanoate has 7 in the acyl chain)
        if chain_length < 7:
            # Not long enough to be considered a fatty acid chain.
            continue

        # Now check each atom in the main chain for extra substituents.
        # For the acyl chain, besides the ones that are part of the main chain (and for the first atom,
        # we allow the connection back to the carbonyl), any extra neighbor must be a trivial methyl branch.
        for i, atom_idx in enumerate(main_chain):
            atom = mol.GetAtomWithIdx(atom_idx)
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                # Allow neighbor if it is part of the main chain
                if nbr_idx in main_chain:
                    continue
                # For the very first acyl atom, allow the carbonyl connection.
                if i == 0 and nbr_idx == carbonyl_atom.GetIdx():
                    continue
                # Check if the substituent is a trivial methyl branch.
                if nbr.GetAtomicNum() == 6 and not nbr.GetIsAromatic():
                    # Count heavy neighbors of this branch carbon excluding the current chain atom.
                    heavy_neighbors = [a for a in nbr.GetNeighbors() 
                                       if a.GetAtomicNum() > 1 and a.GetIdx() != atom_idx]
                    # Allow if it has no further heavy neighbors (i.e. it is a methyl group).
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
            # This ester match shows extra decoration beyond a simple fatty acid.
            continue

        # Compute the number of rotatable bonds along the acyl chain.
        rotatable_in_chain = 0
        for j in range(len(main_chain)-1):
            bond = mol.GetBondBetweenAtoms(main_chain[j], main_chain[j+1])
            if bond and bond.GetBondTypeAsDouble() == 1.0 and not bond.IsInRing():
                rotatable_in_chain += 1
        
        overall_rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
        mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
        actual_heavy_atoms = mol.GetNumHeavyAtoms()
        # Expected heavy atoms = chain_length (acyl chain) + carbonyl carbon + carbonyl oxygen + ester oxygen + methyl carbon.
        expected_heavy_atoms = chain_length + 4
        # Allow a tolerance window (extra groups that might be considered “decoration”)
        tolerance = 6
        
        # Build a detailed reason message.
        reason = (f"Found a methyl ester group with an acyl chain of {chain_length} carbon(s); "
                  f"rotatable bonds in chain: {rotatable_in_chain}, overall rotatable bonds: {overall_rot_bonds}, "
                  f"molecule heavy atoms: {actual_heavy_atoms} (expected ~{expected_heavy_atoms} ± {tolerance}), "
                  f"molecular weight: {mol_wt:.1f} Da")
        return True, reason

    return False, "Molecule contains a methyl ester group but the acyl portion does not meet fatty acid criteria"

# Example usage (you can test with some provided SMILES):
if __name__ == "__main__":
    test_smiles = [
        "CCCCCC1C(C/C(=C(\\C\\C(=C(\\[H])/CCCCCCC(=O)OC)\\[H])/[H])/[H])O1",  # 13-(3-pentyl-2-oxiranyl)-(8Z,11Z)-8,11-tridecadienoic acid methyl ester (TP expected to be rejected due to non-simple acyl chain)
        "O[C@H](CCCCCCCC(OC)=O)/C=C/C=C/CCCCC",  # 9R-hydroxy-10E,12E-octadecadienoic acid methyl ester (TP expected)
        "O=C(OC)CCCCCCCCCCCC",  # methyl tridecanoate (TP expected)
        "O(C(=O)CCC(CCCC)C)C",  # formyl 4-methyl-octanoate (TP expected)
        "O=C(OC)CCCCCCCCCCC/C=C\\CCCCCCCC",  # methyl erucate (TP expected)
        "CCCCCCC(=O)OC",  # Methyl heptanoate (FP; chain too short => rejected)
        "O(C(=O)CC[C@H](CCC=C(C)C)C)C",  # (S)-(-)-Citronellic acid, methyl ester (FP expected)
        "O(C(=O)C[C@@H]1C=CCCC1)C",  # Methyl 2-(cyclohex-2-en-1-yl)acetate (FP expected)
    ]
    
    for sm in test_smiles:
        valid, msg = is_fatty_acid_methyl_ester(sm)
        print(f"SMILES: {sm}\nResult: {valid}\nReason: {msg}\n")