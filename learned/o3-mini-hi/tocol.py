"""
Classifies: CHEBI:39437 tocol
"""
"""
Classifies: tocol
Definition: A chromanol (2,3-dihydro-1-benzopyran) with a chroman-6-ol skeleton that is substituted at the 2‑position 
            by a saturated or triply-unsaturated hydrocarbon chain (three isoprenoid units).

This implementation first identifies a fused bicyclic system in two parts:
  1. A saturated six-membered (dihydro) ring containing one oxygen (the pyran half)
  2. An aromatic six-membered ring (benzene) sharing at least two atoms with the saturated ring.
Then, it checks that the aromatic ring carries a phenolic –OH.
Finally, it identifies a side chain attached to an atom of the saturated ring that is not in the shared aromatic moiety,
ensuring that the side chain is aliphatic (only carbons), sufficiently long, and either saturated (0 double bonds)
or exactly triply unsaturated (3 double bonds).
Note: This heuristic approach is designed to improve over the previous fixed SMARTS method.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdchem import BondType

def is_tocol(smiles: str):
    """
    Determines if the molecule given as a SMILES string is a tocol.
    
    Args:
      smiles (str): SMILES string representing the molecule.
      
    Returns:
      bool: True if the molecule fits the tocol definition, False otherwise.
      str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information from the molecule.
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found in the molecule"
    
    # Find candidate saturated rings (6-membered with exactly one oxygen atom).
    candidate_saturated_rings = []
    for ring in ring_info:
        if len(ring) != 6:
            continue
        atom_types = [mol.GetAtomWithIdx(idx).GetAtomicNum() for idx in ring]
        # exactly one oxygen and rest carbons
        if atom_types.count(8) != 1 or atom_types.count(6) != 5:
            continue
        # Ensure these atoms are not aromatic (the dihydro ring should be non‐aromatic)
        if any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        candidate_saturated_rings.append(set(ring))
    
    if not candidate_saturated_rings:
        return False, "No suitable saturated 6-membered ring containing one oxygen found (dihydro part)"
    
    # Find candidate aromatic rings (benzene) from ring info.
    candidate_aromatic_rings = []
    for ring in ring_info:
        if len(ring) != 6:
            continue
        # Check if all atoms in the ring are aromatic (a benzene ring typically)
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            candidate_aromatic_rings.append(set(ring))
    
    if not candidate_aromatic_rings:
        return False, "No aromatic 6-membered ring (benzene moiety) found"
    
    # Now try to find a fused bicyclic system: a candidate saturated ring fused with an aromatic ring.
    found_core = False
    core_atoms = None
    core_saturated = None
    core_aromatic = None
    core_reason = ""
    
    for sat_ring in candidate_saturated_rings:
        for aro_ring in candidate_aromatic_rings:
            # They must share at least two atoms
            if len(sat_ring.intersection(aro_ring)) >= 2:
                # Build the core as union of the two rings.
                core_atoms = sat_ring.union(aro_ring)
                core_saturated = sat_ring
                core_aromatic = aro_ring
                found_core = True
                break
        if found_core:
            break
            
    if not found_core:
        return False, "No fused bicyclic chromanol core (dihydro ring fused to benzene) found"
    
    # Check that the aromatic (benzene) ring carries a phenolic OH.
    # We scan all atoms in the aromatic ring and look for a neighbor that is an oxygen with at least one hydrogen.
    phenolic_found = False
    for idx in core_aromatic:
        atom = mol.GetAtomWithIdx(idx)
        if not atom.GetIsAromatic():
            continue
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                # Check if this oxygen is bonded to at least one hydrogen.
                # (We look at explicit hydrogens; note that implicit hydrogens are not in the atom list.)
                # We also allow cases where the oxygen belongs to the core itself.
                if nbr.GetTotalNumHs() > 0:
                    phenolic_found = True
                    break
        if phenolic_found:
            break
    if not phenolic_found:
        return False, "No phenolic hydroxyl (-OH) substituent found on the aromatic part of the core"
    
    # Identify the side chain at the 2‑position.
    # Our strategy: in the saturated (dihydro) ring, choose atoms that are NOT shared with the aromatic ring.
    # Then, see if any of these atoms has a neighbor outside the entire core.
    side_chain_root = None
    sat_only = core_saturated.difference(core_aromatic)
    for idx in sat_only:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in core_atoms:
                continue
            if nbr.GetAtomicNum() == 6:  # expecting a carbon
                side_chain_root = nbr
                break
        if side_chain_root is not None:
            # Assume the first valid one is the side chain attachment point.
            break
    if side_chain_root is None:
        return False, "No side chain attached to the non-aromatic portion of the saturated ring (expected at 2‑position)"
    
    # Extract the side chain: perform a breadth-first traversal starting from the side chain root,
    # excluding any atoms that are part of the core.
    side_chain_idxs = set()
    to_visit = [side_chain_root.GetIdx()]
    while to_visit:
        curr = to_visit.pop()
        if curr in side_chain_idxs:
            continue
        side_chain_idxs.add(curr)
        for nbr in mol.GetAtomWithIdx(curr).GetNeighbors():
            if nbr.GetIdx() in core_atoms:
                continue
            if nbr.GetIdx() not in side_chain_idxs:
                to_visit.append(nbr.GetIdx())
    
    if not side_chain_idxs:
        return False, "Failed to extract the 2‑position side chain"
    
    # Check that the side chain contains only carbon atoms.
    for idx in side_chain_idxs:
        if mol.GetAtomWithIdx(idx).GetAtomicNum() != 6:
            return False, f"Side chain contains a non-carbon atom ({mol.GetAtomWithIdx(idx).GetSymbol()} at index {idx})"
    
    n_chain_carbons = len(side_chain_idxs)
    if n_chain_carbons < 10:
        return False, f"Side chain is too short (found {n_chain_carbons} carbons; expected at least 10)"
    
    # Count double bonds within the side chain.
    double_bond_count = 0
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        if a1 in side_chain_idxs and a2 in side_chain_idxs:
            if bond.GetBondType() == BondType.DOUBLE:
                double_bond_count += 1

    if double_bond_count not in (0, 3):
        return False, f"Side chain has {double_bond_count} double bond(s); expected 0 (fully saturated) or 3 (triply unsaturated)"
    
    return True, "Molecule contains a chromanol core (fused dihydropyran and benzene with phenolic -OH) with a proper 2‑position side chain"

# For testing (you can remove or modify these lines as needed):
if __name__ == "__main__":
    test_smiles_list = [
        "CC(C)CCC[C@H](C)CCC[C@@H](C)CCC[C@@]1(C)CCC2=C(C)C(O)=C(C)C(C)=C2O1",  # (S,R,S)-alpha-tocopherol (expected to be tocol)
        "Cl.NC[C@@H]1O[C@@H](Cc2c(O)c(O)ccc12)C12CC3CC(CC(C3)C1)C2",  # False positive example reported earlier
        "O1C(CCCCCCCCCCCCC)CC2=C(C1=O)C(O)=CC=C2",  # Another false positive candidate
    ]
    for s in test_smiles_list:
        flag, reason = is_tocol(s)
        print("SMILES:", s)
        print("Result:", flag, "|", reason)
        print("------------------------------------------------")