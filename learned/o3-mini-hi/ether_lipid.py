"""
Classifies: CHEBI:64611 ether lipid
"""
"""
Classifies: Ether Lipid
Definition: A lipid similar in structure to a glycerolipid but in which one or more 
of the carbon atoms on glycerol is bonded to an alkyl chain via an ether linkage (C–O–C)
rather than the usual ester linkage.
Heuristic criteria used in this version:
  (1) The molecule must have a molecular weight typical for a lipid (>= 300 Da).
  (2) Search for an ether oxygen (an oxygen atom with exactly two carbon neighbors)
      that is not part of a carbonyl group.
  (3) For each such oxygen, consider both attached carbons:
      a. One (the “polar candidate”) must be embedded in a small glycerol‐like fragment.
         In our heuristic, we require that aside from the ether oxygen it is attached to at least one 
         hydroxyl group itself and that one of its other carbon neighbors also carries an –OH.
      b. The other (the “chain candidate”) must be acyclic and lead to a continuous alkyl chain 
         (through single C–C bonds) of at least 8 carbons.
If any such ether bond is found, the molecule is classified as an ether lipid.
Note that these heuristics are approximate.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as an ether lipid, else False.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Step 1: Ensure that the molecular weight is above the typical lower bound for lipids.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 300:
        return False, "Molecular weight too low for a lipid"

    # Helper: Check if a carbon atom (candidate for linkage) is part of a carbonyl group.
    def is_carbonyl(carbon):
        for bond in carbon.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                other = bond.GetOtherAtom(carbon)
                if other.GetAtomicNum() == 8:
                    return True
        return False

    # Helper: Recursively compute the length of the longest contiguous, acyclic carbon chain 
    # starting from the given atom. Only count single bonds and skip atoms in rings.
    def longest_chain(atom, visited):
        if atom.IsInRing():
            return 0
        max_length = 1  # count self
        visited.add(atom.GetIdx())
        for bond in atom.GetBonds():
            if bond.GetBondType() != Chem.BondType.SINGLE:
                continue
            nbr = bond.GetOtherAtom(atom)
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited and not nbr.IsInRing():
                branch_length = 1 + longest_chain(nbr, visited.copy())
                if branch_length > max_length:
                    max_length = branch_length
        return max_length

    # Helper: Check if an oxygen atom is a hydroxyl (i.e. –OH) by verifying it has at least one H.
    def is_hydroxyl(ox):
        # Using GetTotalNumHs() gives the total (implicit + explicit) number of hydrogens.
        return ox.GetAtomicNum() == 8 and ox.GetTotalNumHs() >= 1

    # New helper: Evaluate whether a candidate polar carbon (attached to the ether oxygen)
    # is present in a glycerol-like fragment.
    # Heuristically, require that (a) it has at least one hydroxyl neighbor (other than the ether oxygen)
    # and (b) at least one of its other carbon neighbors carries a hydroxyl group.
    def in_glycerol_fragment(polar_candidate, ether_oxygen):
        hydroxyl_count = 0
        # Check direct neighbors of the polar_candidate (skip the ether oxygen)
        for nbr in polar_candidate.GetNeighbors():
            if nbr.GetIdx() == ether_oxygen.GetIdx():
                continue
            if is_hydroxyl(nbr):
                hydroxyl_count += 1
        # Now check if any carbon neighbor has an –OH substituent.
        for nbr in polar_candidate.GetNeighbors():
            if nbr.GetIdx() == ether_oxygen.GetIdx() or nbr.GetAtomicNum() != 6:
                continue
            for nbr2 in nbr.GetNeighbors():
                if nbr2.GetIdx() in (polar_candidate.GetIdx(), ether_oxygen.GetIdx()):
                    continue
                if is_hydroxyl(nbr2):
                    hydroxyl_count += 1
                    break
        return hydroxyl_count >= 2

    # Loop over all oxygen atoms to identify candidate ether oxygens.
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
        # Skip if either carbon is involved in a carbonyl (to avoid ester linkages).
        if is_carbonyl(c1) or is_carbonyl(c2):
            continue

        # Try both orderings:
        for polar_candidate, chain_candidate in [(c1, c2), (c2, c1)]:
            # For the polar side, require at least one hydroxyl neighbor (other than the ether oxygen).
            extra_oh = False
            for nbr in polar_candidate.GetNeighbors():
                if nbr.GetIdx() == atom.GetIdx():
                    continue
                if is_hydroxyl(nbr):
                    extra_oh = True
                    break
            if not extra_oh:
                continue  # This candidate does not provide enough polarity.
            
            # Also verify that the polar candidate is in a glycerol-like fragment.
            if not in_glycerol_fragment(polar_candidate, atom):
                continue

            # The chain candidate must not be in a ring (to ensure a contiguous chain).
            if chain_candidate.IsInRing():
                continue

            # Now evaluate the length of the alkyl chain stemming from the chain candidate.
            chain_length = longest_chain(chain_candidate, set())
            if chain_length >= 8:
                return (True, 
                        "Molecule contains a glycerol-like headgroup with at least two hydroxyls and an ether linkage "
                        "to a long alkyl chain (chain length {} >= 8)".format(chain_length))
    
    return False, "No suitable ether linkage with a proper glycerol-like headgroup and a long alkyl chain found"


# Example usage (testing a few SMILES strings):
if __name__ == "__main__":
    test_smiles = [
        "P(OC[C@H](OC(=O)CCCCCCCCC/C=C\\CCCCCCCCCC)COCCCCCCCCCCCCCCCCCC)(O)(O)=O",  # candidate (PA(O-18:0/22:1(11Z)))
        "C[C@H]1CCC[C@@H](C)CCC[C@@H](C)CC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCOC(CO)COCC[C@@H](C)CCC1",  # caldarchaeol candidate
        "CCCCCCCC\\C=C/CCCCCCCCOC[C@@H](O)COP([O-])(=O)OCC[N+](C)(C)C",  # candidate 1-oleyl-sn-glycero-3-phosphocholine
        "C([C@@](CO/C=C\\CCCCCCCCCCCCCCCC)(O)[H])OP([O-])(=O)OCC[N+](C)(C)C"  # non-ether lipid example
    ]
    for s in test_smiles:
        res, reason = is_ether_lipid(s)
        print("SMILES:", s)
        print("Result:", res)
        print("Reason:", reason)
        print()