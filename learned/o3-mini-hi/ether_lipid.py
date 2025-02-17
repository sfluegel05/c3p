"""
Classifies: CHEBI:64611 ether lipid
"""
"""
Classifies: Ether Lipid
Definition: A lipid similar in structure to a glycerolipid but in which one or more 
of the carbon atoms on glycerol is bonded to an alkyl chain via an ether linkage (C–O–C)
rather than the usual ester linkage.
Heuristic criteria used:
  (1) The molecule must have a molecular weight typical for a lipid (>= 300 Da).
  (2) Look for a candidate ether oxygen (an O with exactly two carbon neighbors) that is not part of a carbonyl.
  (3) For each candidate, treat one of the attached carbons as the “polar candidate” and the other as the “chain candidate”.
      • The polar candidate must have at least one extra oxygen neighbor AND be in a glycerol‐like local environment,
        meaning that at least one of its carbon neighbors (other than the one from the ether bond) also has an oxygen neighbor.
      • The chain candidate must be acyclic (not in a ring) and attached to a contiguous alkyl chain of at least 8 carbons.
If any such ether bond is found, the molecule is classified as an ether lipid.
Note that this implementation is heuristic and was refined in response to many false positives 
(from molecules that only incidentally contained a glycerol‐like substructure and an ether linkage).
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is classified as an ether lipid, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1: ensure the molecular weight is above 300 Da.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 300:
        return False, "Molecular weight too low for a lipid"
    
    # Helper: check if a carbon atom is part of a carbonyl group (has a double bond to oxygen)
    def is_carbonyl(carbon):
        for bond in carbon.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                other = bond.GetOtherAtom(carbon)
                if other.GetAtomicNum() == 8:
                    return True
        return False

    # Helper: compute longest contiguous acyclic chain (only through carbons not in a ring)
    def longest_chain(atom, visited):
        if atom.IsInRing():
            return 0
        max_length = 1  # count self
        visited.add(atom.GetIdx())
        for bond in atom.GetBonds():
            if bond.GetBondType() != Chem.BondType.SINGLE:
                continue
            nbr = bond.GetOtherAtom(atom)
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                if nbr.IsInRing():
                    continue
                branch_length = 1 + longest_chain(nbr, visited.copy())
                if branch_length > max_length:
                    max_length = branch_length
        return max_length

    # Extra glycerol-context check: given a candidate polar carbon (from the ether bond),
    # we expect that in a glycerol backbone at least one of its carbon neighbors (besides the bonded oxygen)
    # is also substituted by an oxygen (i.e. part of a 3‐carbon chain with multiple oxygens).
    def in_glycerol_context(polar_candidate, ether_oxygen):
        # Look at neighbors of polar_candidate (skip the ether oxygen)
        for nbr in polar_candidate.GetNeighbors():
            if nbr.GetIdx() == ether_oxygen.GetIdx():
                continue
            if nbr.GetAtomicNum() != 6:
                continue
            # Check if this neighbor has any oxygen (other than polar_candidate) attached.
            for nbr2 in nbr.GetNeighbors():
                if nbr2.GetIdx() == polar_candidate.GetIdx():
                    continue
                if nbr2.GetAtomicNum() == 8:
                    return True
        return False

    # Now go through every oxygen in the molecule...
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:
            continue
        # For an ether oxygen we expect exactly two neighbors.
        if atom.GetDegree() != 2:
            continue
        nbrs = atom.GetNeighbors()
        # Both neighbors must be carbons.
        if not all(nbr.GetAtomicNum() == 6 for nbr in nbrs):
            continue
        c1, c2 = nbrs[0], nbrs[1]
        # Skip if either carbon is part of a carbonyl (to avoid ester linkages).
        if is_carbonyl(c1) or is_carbonyl(c2):
            continue
        
        # Try both orderings: one as polar candidate (glycerol side) and the other as chain candidate.
        for polar_candidate, chain_candidate in [(c1, c2), (c2, c1)]:
            # Check that polar candidate has at least one extra oxygen neighbor besides the ether oxygen.
            extra_oxygen = False
            for nbr in polar_candidate.GetNeighbors():
                if nbr.GetIdx() == atom.GetIdx():
                    continue
                if nbr.GetAtomicNum() == 8:
                    extra_oxygen = True
                    break
            if not extra_oxygen:
                continue  # Not enough polarity on this side.
            
            # Additional glycerol context requirement.
            if not in_glycerol_context(polar_candidate, atom):
                continue  # This polar candidate does not appear to be part of a glycerol-like group.
                
            # The chain candidate must not be part of a ring.
            if chain_candidate.IsInRing():
                continue
            
            # Search for a contiguous acyclic alkyl chain attached to the chain candidate.
            longest = 0
            for nbr in chain_candidate.GetNeighbors():
                if nbr.GetIdx() == atom.GetIdx():
                    continue
                if nbr.GetAtomicNum() == 6:
                    chain_length = longest_chain(nbr, set())
                    if chain_length > longest:
                        longest = chain_length
            # Account for the chain candidate itself if needed.
            longest = max(longest, 1)
            if longest >= 8:
                return (True, 
                        "Molecule contains a glycerol-like headgroup (satisfying extra oxygen and glycerol context) "
                        "and an ether linkage to a long alkyl chain (chain length {} >=8)".format(longest))
    
    return False, "No suitable ether linkage with a glycerol-like headgroup and a long alkyl chain found"

# Example usage (testing with a small set of SMILES from provided examples):
if __name__ == "__main__":
    test_smiles = [
        "P(OC[C@H](OC(=O)CCCCCCCCC/C=C\\CCCCCCCCCC)COCCCCCCCCCCCCCCCCCC)(O)(O)=O",  # ether lipid candidate
        "C[C@H]1CCC[C@@H](C)CCC[C@@H](C)CC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCOC(CO)COCC[C@@H](C)CCC1",  # caldarchaeol - expected to be ether lipid if glycerol detected
        "CCCCCCCC\\C=C/CCCCCCCCOC[C@@H](O)COP([O-])(=O)OCC[N+](C)(C)C",  # 1-oleyl-sn-glycero-3-phosphocholine candidate
        "C([C@@](CO/C=C\\CCCCCCCCCCCCCCCC)(O)[H])OP([O-])(=O)OCC[N+](C)(C)C"  # non-ether lipid example
    ]
    for s in test_smiles:
        res, reason = is_ether_lipid(s)
        print("SMILES:", s)
        print("Result:", res)
        print("Reason:", reason)
        print()