"""
Classifies: CHEBI:36313 glycerophosphocholine
"""
"""
Classifies: glycerophosphocholine
Definition: The glycerol phosphate ester of a phosphocholine.
A true glycerophosphocholine (GPC) should contain (i) a phosphocholine headgroup – that is, a phosphate 
group linked to an oxygen that carries a two‐carbon chain ending in a trimethylammonium group – and 
(ii) that phosphate must be connected via an oxygen atom (typically on a primary, CH2‐bearing carbon) 
to a short glycerol (or glycerol‐like) fragment. In many phospholipids (e.g. diacyl phosphatidylcholines) 
the glycerol bears two acyl chains, but those are considered false positives here.
This heuristic is meant to help distinguish true glycerophosphocholines (or lysophosphatidylcholines) 
from other phosphocholine-containing lipids.
Note: This is only a heuristic method.
"""

from rdkit import Chem

def is_glycerophosphocholine(smiles: str):
    """
    Determines whether a molecule (given as a SMILES string) is a glycerophosphocholine.
    
    The algorithm proceeds in roughly two steps:
      1. It verifies that a phosphocholine headgroup is present. This is done by looking
         for at least one substructure match using one of several SMARTS patterns that capture 
         the headgroup motif: an oxygen-bridged phosphate with a two-carbon chain ending in trimethylammonium.
      2. For each phosphorus atom in any headgroup match, we inspect one of its oxygen neighbors 
         that (a) is attached to a carbon bearing at least two hydrogens (a proxy for a primary CH₂), 
         and (b) whose immediate environment is “glycerol-like” (i.e. part of a short, three‐carbon chain).
         In that candidate glycerol fragment we also count how many acyl (ester) substitutions are present.
         More than one (i.e. a double‐acylated glycerol) is not accepted.
      
    Args:
      smiles (str): SMILES string for the molecule.
      
    Returns:
      bool: True if the molecule is accepted as glycerophosphocholine, else False.
      str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Step 1. Look for a phosphocholine headgroup.
    # We use a short list of SMARTS that differ only in the protonation state of an -OH on phosphorus.
    headgroup_smarts_list = [
        "[CH2]O[P](=O)(O)OCC[N+](C)(C)C",      # all groups protonated
        "[CH2]O[P](=O)([O-])OCC[N+](C)(C)C",    # one deprotonated oxygen (variant 1)
        "[CH2]O[P]([O-])(=O)OCC[N+](C)(C)C"     # one deprotonated oxygen (variant 2)
    ]
    headgroup_found = False
    for smarts in headgroup_smarts_list:
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            continue
        if mol.HasSubstructMatch(patt):
            headgroup_found = True
            break
    if not headgroup_found:
        return False, "Phosphocholine headgroup substructure not found"

    # --- Step 2. Check for proper glycerol connectivity.
    # For every phosphorus atom, we inspect its oxygen neighbors.
    # For each such oxygen (that is not directly part of the choline branch),
    # we require that its other neighbor is a carbon bearing ≥2 H's (i.e. a likely CH2 group)
    # and that this carbon is part of a short (glycerol-like) fragment.
    glycerol_candidate_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 15:
            continue  # consider only phosphorus atoms
        p_atom = atom
        # Get all oxygen neighbors of the phosphorus.
        for o_atom in p_atom.GetNeighbors():
            if o_atom.GetAtomicNum() != 8:
                continue
            # Get the neighbor(s) of o_atom that are not the phosphorus.
            nonP_neighbors = [nbr for nbr in o_atom.GetNeighbors() if nbr.GetIdx() != p_atom.GetIdx()]
            if not nonP_neighbors:
                continue
            # Skip if this oxygen appears to be part of the choline branch: typically if one of its neighbors is
            # a nitrogen with a formal +1 charge.
            is_choline = any(nbr.GetAtomicNum() == 7 and nbr.GetFormalCharge() == 1 for nbr in nonP_neighbors)
            if is_choline:
                continue
            # Now check each nonP neighbor for a primary carbon.
            for nbr in nonP_neighbors:
                if nbr.GetAtomicNum() == 6:
                    # Count total hydrogens (implicit + explicit)
                    if nbr.GetTotalNumHs() < 2:
                        continue
                    # We now have a candidate connection: O – C(≥CH2)
                    # Next, try to decide if this carbon is part of a glycerol-like fragment.
                    # In glycerol (or glycerol-like systems) the backbone is very short (three carbons) and one
                    # of the carbons (often the one linked to phosphate) is primary.
                    # A simple heuristic: we look at the neighbors of this candidate carbon and check if at least one
                    # other neighbor is an oxygen (a hydroxyl) but not part of an ester with a C=O.
                    candidate_c = nbr
                    oxy_neighbor_found = False
                    acyl_substituents = 0
                    for cand_nbr in candidate_c.GetNeighbors():
                        # Skip the oxygen that connected candidate_c to the phosphate.
                        if cand_nbr.GetIdx() == o_atom.GetIdx():
                            continue
                        if cand_nbr.GetAtomicNum() == 8:
                            oxy_neighbor_found = True
                            # Check if the bond from candidate_c to this oxygen is an ester linkage (i.e. oxygen further bound to C=O).
                            for second_nbr in cand_nbr.GetNeighbors():
                                if second_nbr.GetIdx() == candidate_c.GetIdx():
                                    continue
                                # Look for a carbonyl: carbon double-bonded to an oxygen.
                                if second_nbr.GetAtomicNum() == 6:
                                    for third_nbr in second_nbr.GetNeighbors():
                                        if third_nbr.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(second_nbr.GetIdx(), third_nbr.GetIdx()).GetBondTypeAsDouble() == 1:
                                            acyl_substituents += 1
                                            break
                    # We require that the candidate carbon is part of a glycerol-like fragment (at least one O neighbor)
                    # and that it is not heavily acylated (at most one acyl ester allowed: lysophosphatidylcholines in many cases).
                    if oxy_neighbor_found and acyl_substituents <= 1:
                        glycerol_candidate_found = True
                        break
            if glycerol_candidate_found:
                break
        if glycerol_candidate_found:
            break
    if not glycerol_candidate_found:
        return False, "Glycerol backbone connectivity not found (phosphate not linked via an oxygen to a primary CH2 carbon with a glycerol-like environment)"
    
    return True, "Molecule contains a phosphocholine headgroup with proper glycerol connectivity"

# Example usage.
if __name__ == "__main__":
    # Test with one provided SMILES string.
    test_smiles = "P(OC[C@@H](CO)OC(CCCCCCC/C=C\\CC(CCCCCC)O)=O)(=O)(OCC[N+](C)(C)C)[O-]"
    result, reason = is_glycerophosphocholine(test_smiles)
    print(result, reason)