"""
Classifies: CHEBI:64611 ether lipid
"""
"""
Classifies: Ether Lipid
Definition: A lipid similar in structure to a glycerolipid but in which one or more 
of the carbon atoms on glycerol is bonded to an alkyl chain via an ether linkage (C–O–C)
rather than the usual ester linkage.
Heuristic criteria for this version:
  (1) The molecule must have a molecular weight typical for lipids (>= 300 Da).
  (2) Loop over oxygen atoms with exactly two carbon neighbors.
  (3) Exclude ethers where either attached carbon is involved in a carbonyl (i.e. part of an ester).
  (4) For each candidate oxygen, try both orders:
      - One attached carbon is considered the "polar candidate" (glycerol‐like headgroup).
         We now require this carbon to have at least one neighboring oxygen or phosphorus
         (other than the ether oxygen).
      - The other attached carbon is the "chain candidate" and must be acyclic and lead
         to a linear alkyl chain of at least 8 carbons in length.
If any candidate pair satisfies these conditions, the molecule is classified as an ether lipid.
Note that these heuristics are approximate.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an ether lipid, else False.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Ensure that the molecular weight is above a typical lower bound for lipids.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 300:
        return False, f"Molecular weight too low for a lipid ({mw:.1f} < 300 Da)"
    
    # Helper: Check if a carbon atom is part of a carbonyl group (C=O)
    def is_carbonyl(carbon):
        for bond in carbon.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                other = bond.GetOtherAtom(carbon)
                if other.GetAtomicNum() == 8:
                    return True
        return False

    # Helper: For a given starting carbon atom, recursively find the longest contiguous acyclic carbon chain.
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

    # New helper: Check if the polar candidate (attached carbon) is in a glycerol-like environment.
    # Instead of requiring two hydroxyl groups, we now check if it is substituted with at least one
    # oxygen (or phosphorus) that is not the ether oxygen.
    def in_glycerol_like(polar_candidate, ether_oxygen):
        for nbr in polar_candidate.GetNeighbors():
            if nbr.GetIdx() == ether_oxygen.GetIdx():
                continue
            # Accept if neighbor is oxygen (atomic num 8) or phosphorus (atomic num 15)
            if nbr.GetAtomicNum() in (8, 15):
                return True
        return False

    # Loop over all oxygen atoms in the molecule to pick candidate ether linkages.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:
            continue
        # For an ether oxygen, we expect exactly two neighbors.
        if atom.GetDegree() != 2:
            continue
        nbrs = atom.GetNeighbors()
        # Both neighbors must be carbons.
        if not all(nbr.GetAtomicNum() == 6 for nbr in nbrs):
            continue

        c1, c2 = nbrs[0], nbrs[1]
        # Exclude this oxygen if either attached carbon is in a carbonyl (to avoid misidentifying ester bonds)
        if is_carbonyl(c1) or is_carbonyl(c2):
            continue

        # Try both orderings: assign one as polar candidate and the other as chain candidate.
        for polar_candidate, chain_candidate in [(c1, c2), (c2, c1)]:
            # Requirement for the polar candidate: should have at least one oxygen/phosphorus neighbor (aside from the ether oxygen)
            if not in_glycerol_like(polar_candidate, atom):
                continue

            # Requirement for the chain candidate: must not be in a ring
            if chain_candidate.IsInRing():
                continue

            # Measure the longest chain starting from the chain candidate.
            chain_length = longest_chain(chain_candidate, set())
            if chain_length >= 8:
                return (True, 
                        "Molecule contains at least one ether linkage: a polar (glycerol-like) side with an oxygen/phosphorus "
                        "substituent and a long alkyl chain (chain length {} >= 8)".format(chain_length))
    
    return False, "No suitable ether linkage with a proper glycerol-like headgroup and a long alkyl chain found"


# Example usage (testing a few SMILES strings):
if __name__ == "__main__":
    test_smiles = [
        "P(OC[C@H](OC(=O)CCCCCCCCC/C=C\\CCCCCCCCCC)COCCCCCCCCCCCCCCCCCC)(O)(O)=O",  # PA(O-18:0/22:1(11Z))
        "C[C@H]1CCC[C@@H](C)CCC[C@@H](C)CC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCOC(CO)COCC[C@@H](C)CCC[C@H](C)CCC[C@@H](C)CCC[C@@H](C)CC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCOC(CO)COCC[C@@H](C)CCC1",  # caldarchaeol candidate
        "CCCCCCCC\\C=C/CCCCCCCCOC[C@@H](O)COP([O-])(=O)OCC[N+](C)(C)C",  # 1-oleyl-sn-glycero-3-phosphocholine candidate
        "C([C@@](CO/C=C\\CCCCCCCCCCCCCCCC)(O)[H])OP([O-])(=O)OCC[N+](C)(C)C"  # Non-ether lipid example
    ]
    for s in test_smiles:
        res, reason = is_ether_lipid(s)
        print("SMILES:", s)
        print("Result:", res)
        print("Reason:", reason)
        print()