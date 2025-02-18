"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
"""
Classifies: Iridoid monoterpenoid
Definition: Monoterpenoids biosynthesized from isoprene units. Typically they contain
a cyclopentane ring fused to a six‐membered oxygen heterocycle (the “iridoid core”).
Cleavage of a bond in the cyclopentane ring gives rise to secoiridoids.
This program uses two strategies:
  (1) Look for a fused bicyclic system where a 5‐membered ring and a 6‐membered ring share
      at least 2 atoms. Then merge the rings (union) – if the total heavy‐atom count is 9
      with 8 carbons and 1 oxygen (the “ideal” iridoid core), we classify as iridoid.
  (2) If not found, then try to “rescue” secoiridoids by identifying a cyclopentane ring (5 carbons)
      that carries at least two exocyclic carbonyl groups (C=O on a neighbor atom).
      
Note: This heuristic may still mis‐classify some molecules.
"""

from rdkit import Chem

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines if a molecule is an iridoid monoterpenoid (or its secoiridoid derivative)
    based on its SMILES string.
    
    It works in two stages:
      Stage 1: Look for an ideal fused bicyclic core composed of a 5-membered ring (cyclopentane)
               fused with a 6-membered ring that contains one oxygen. When the atoms in the two rings
               are merged (shared atoms counted only once) the core should have 9 heavy atoms:
               8 carbons and 1 oxygen.
      Stage 2: If no fused bicyclic core is found, search for a cyclopentane ring (5-membered ring
               composed solely of carbons) that has at least two exocyclic carbonyl (C=O) groups.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an iridoid monoterpenoid (or secoiridoid),
              False otherwise.
        str: Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found in the molecule"
    
    # Stage 1: Look for fused bicyclic core with a 5-membered and a 6-membered ring sharing >=2 atoms.
    for i in range(len(ring_info)):
        for j in range(i+1, len(ring_info)):
            ring1 = set(ring_info[i])
            ring2 = set(ring_info[j])
            # Check for one ring of size 5 and the other of size 6.
            if (len(ring1) == 5 and len(ring2) == 6) or (len(ring1) == 6 and len(ring2) == 5):
                shared_atoms = ring1.intersection(ring2)
                if len(shared_atoms) >= 2:
                    # Get the union of the two rings.
                    core = ring1.union(ring2)
                    # We expect an ideal iridoid core to have 9 heavy atoms.
                    if len(core) == 9:
                        # Count carbons and oxygens in the fused core.
                        nC = 0
                        nO = 0
                        for idx in core:
                            atom = mol.GetAtomWithIdx(idx)
                            if atom.GetAtomicNum() == 6:
                                nC += 1
                            elif atom.GetAtomicNum() == 8:
                                nO += 1
                        if nC == 8 and nO == 1:
                            reason = ("Found fused bicyclic core: a 5-membered ring and a 6-membered oxygen heterocycle "
                                      "merging to 9 atoms (8 carbons, 1 oxygen) with {} shared atoms.".format(len(shared_atoms)))
                            return True, reason
                            
    # Stage 2: Rescue approach for secoiridoids:
    # Look for a cyclopentane ring (5-membered) that is composed entirely of carbons.
    # Then count how many exocyclic carbonyl groups (C=O) are attached to atoms of the ring.
    for ring in ring_info:
        if len(ring) == 5:
            is_cyclopentane = True
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() != 6:
                    is_cyclopentane = False
                    break
            if not is_cyclopentane:
                continue
            carbonyl_count = 0
            # For each atom in the cyclopentane ring, examine neighbors not in the ring.
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                for nbr in atom.GetNeighbors():
                    nbr_idx = nbr.GetIdx()
                    if nbr_idx in ring:
                        continue
                    bond = mol.GetBondBetweenAtoms(idx, nbr_idx)
                    # Check if this neighbor carbon is part of a carbonyl group.
                    # We require the bond to be a double bond and the neighbor to be carbon.
                    if bond is not None and bond.GetBondTypeAsDouble() == 2 and nbr.GetAtomicNum() == 6:
                        # Now check if the neighbor carbon is double-bonded to an oxygen.
                        for nbr2 in nbr.GetNeighbors():
                            if nbr2.GetIdx() == idx:
                                continue
                            bond2 = mol.GetBondBetweenAtoms(nbr.GetIdx(), nbr2.GetIdx())
                            if bond2 is not None and bond2.GetBondTypeAsDouble() == 2 and nbr2.GetAtomicNum() == 8:
                                carbonyl_count += 1
                                break  # Count once per exocyclic carbonyl
            if carbonyl_count >= 2:
                reason = ("Found cyclopentane ring (all carbon atoms) with {} exocyclic carbonyl substituents "
                          "suggestive of a secoiridoid skeleton.".format(carbonyl_count))
                return True, reason

    return False, "No fused iridoid or secoiridoid core detected based on our heuristics."

# Example usage (for testing purposes):
if __name__ == "__main__":
    # Test examples: (true positives and negatives from the provided outcomes)
    examples = [
        # True positive: fused core detected
        ("[C@]12([C@@]([C@H](CC1)C)(C(OC=C2C)=O)[H])[H]", "cis-trans-nepetalactone"),
        # An example secoiridoid: iridodial (one of the missed ones before)
        ("CC1CCC(C(C)C=O)C1C=O", "iridodial"),
        # A false positive example from the outcomes:
        ("C1[C@@H](O[C@@H]([C@H]2[C@@H]1C3=C(O2)C=CC(=C3)NS(=O)(=O)C4=CC=CC=C4)CO)CC(=O)NCC(F)(F)F", "false positive example"),
        # A false negative from prior run: oleuropein aglycone
        ("O[C@H]1/C(/[C@](C(=CO1)C(=O)OC)([H])CC(=O)OCCC=2C=CC(=C(C2)O)O)=C/C", "oleuropein aglycone"),
    ]
    
    for smi, name in examples:
        flag, msg = is_iridoid_monoterpenoid(smi)
        print(f"{name}: {flag} => {msg}")