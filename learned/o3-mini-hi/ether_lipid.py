"""
Classifies: CHEBI:64611 ether lipid
"""
"""
Classifies: Ether Lipid
Definition: A lipid similar in structure to a glycerolipid but in which one or more of 
the carbon atoms on glycerol is bonded to an alkyl chain via an ether linkage (C-O-C) 
rather than the usual ester linkage.
Heuristic criteria used:
  (1) The molecule must have a molecular weight typical for a lipid (>=300 Da).
  (2) There must be at least one ether linkage in which the central oxygen atom is bonded
      exactly to two carbons.
  (3) One side of that ether (the “polar” side) must carry at least one oxygen neighbor 
      other than the ether oxygen (a rough sign of a glycerol-like headgroup).
  (4) The other side (the “chain” side) must be attached to a long contiguous alkyl chain.
      In our search we only follow carbons that are not in any ring (i.e. acyclic chain) 
      and count the number of connected carbons (including the starting one). We require 
      a chain length of at least 8.
  (5) We also skip ether bonds if either carbon is “carbonyl‐like” (i.e. is double‐bonded 
      to an oxygen) to avoid picking up ester linkages.
This implementation is heuristic and may be further refined.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid based on its SMILES string.
    The approach is as follows:
      - Ensure the molecular weight is above 300 Da.
      - For every oxygen atom that has exactly two carbon neighbors (candidate ether oxygen):
            * Check that neither attached carbon is part of a carbonyl.
            * For each ordering of the two attached carbons, designate one as the
              candidate polar (glycerol-like) carbon and the other as the candidate chain carbon.
            * The polar candidate must have at least one additional oxygen neighbor (besides the ether oxygen).
            * The candidate chain side must not be within a ring and must be attached to
              a contiguous chain of (at least 8) carbon atoms.
      - If any such ether bond is found, return True.
      
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as an ether lipid, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Molecules that are too small are unlikely lipids.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 300:
        return False, "Molecular weight too low for a lipid"
    
    # Helper: check if a carbon atom is part of a carbonyl group
    def is_carbonyl(carbon):
        for bond in carbon.GetBonds():
            # Look for a double bond to oxygen
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                other = bond.GetOtherAtom(carbon)
                if other.GetAtomicNum() == 8:
                    return True
        return False

    # Helper: recursively compute the longest contiguous acyclic chain (counting carbon atoms)
    # starting from a given carbon. We only walk through carbon atoms that are not in any ring.
    def longest_chain(atom, visited):
        # If this carbon is in a ring, we do not count it as part of a linear alkyl chain.
        if atom.IsInRing():
            return 0
        max_length = 1  # count self
        visited.add(atom.GetIdx())
        for bond in atom.GetBonds():
            # Only consider single bonds to carbons
            if bond.GetBondType() != Chem.BondType.SINGLE:
                continue
            nbr = bond.GetOtherAtom(atom)
            if nbr.GetAtomicNum() == 6 and (nbr.GetIdx() not in visited):
                # Only follow neighbor if it is not in a ring.
                if nbr.IsInRing():
                    continue
                branch_length = 1 + longest_chain(nbr, visited.copy())
                if branch_length > max_length:
                    max_length = branch_length
        return max_length

    # Now, go through all oxygen atoms to see if any qualifies as the ether oxygen.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:
            continue
        # For an ether oxygen we expect exactly two neighbors.
        if atom.GetDegree() != 2:
            continue
        nbrs = atom.GetNeighbors()
        if not all(nbr.GetAtomicNum() == 6 for nbr in nbrs):
            continue
        c1, c2 = nbrs[0], nbrs[1]
        # Exclude if either carbon is part of a carbonyl group.
        if is_carbonyl(c1) or is_carbonyl(c2):
            continue
        
        # For each ordering, designate one carbon as the "polar" candidate, the other as "chain" candidate.
        for polar_candidate, chain_candidate in [(c1, c2), (c2, c1)]:
            # For the polar candidate, check for an extra oxygen neighbor (besides our ether oxygen).
            extra_oxygen = False
            for nbr in polar_candidate.GetNeighbors():
                if nbr.GetIdx() == atom.GetIdx():
                    continue
                if nbr.GetAtomicNum() == 8:
                    extra_oxygen = True
                    break
            if not extra_oxygen:
                continue  # not polar enough
            
            # For the chain candidate, require that it is not part of a ring.
            if chain_candidate.IsInRing():
                continue
            
            # Search for a contiguous acyclic alkyl chain on the chain candidate.
            longest = 0
            # Look at each neighbor of chain_candidate (except the ether oxygen).
            for nbr in chain_candidate.GetNeighbors():
                if nbr.GetIdx() == atom.GetIdx():
                    continue
                if nbr.GetAtomicNum() == 6:
                    # Calculate chain length starting from this neighbor.
                    chain_length = longest_chain(nbr, set())
                    if chain_length > longest:
                        longest = chain_length
            # Also count the chain candidate itself.
            longest = max(longest, 1)
            if longest >= 8:
                return (True, 
                        "Molecule contains a glycerol-like moiety (carbon with extra oxygen) "
                        "and an ether linkage to a long alkyl chain (chain length {} >=8)".format(longest))
    
    return False, "No suitable ether linkage with a glycerol-like headgroup and a long alkyl chain found"

# Example usage:
if __name__ == "__main__":
    # A small set of test SMILES from the provided examples.
    test_smiles = [
        "P(OC[C@H](OC(=O)CCCCCCCCC/C=C\\CCCCCCCCCC)COCCCCCCCCCCCCCCCCCC)(O)(O)=O",  # expected positive
        "C[C@H]1CCC[C@@H](C)CCC[C@@H](C)CC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCOC(CO)COCC[C@@H](C)CCC1",  # expected missed previously
        "CCCCCCCC\\C=C/CCCCCCCCOC[C@@H](O)COP([O-])(=O)OCC[N+](C)(C)C",  # expected positive ether lipid
        "C([C@@](CO/C=C\\CCCCCCCCCCCCCCCC)(O)[H])OP([O-])(=O)OCC[N+](C)(C)C"  # expected negative
    ]
    for s in test_smiles:
        res, reason = is_ether_lipid(s)
        print("SMILES:", s)
        print("Result:", res)
        print("Reason:", reason)
        print()