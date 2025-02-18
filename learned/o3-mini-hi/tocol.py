"""
Classifies: CHEBI:39437 tocol
"""
"""
Classifies: tocol
Definition: A chromanol (2,3-dihydro-1-benzopyran) with a chroman-6-ol skeleton that is substituted at the 2‑position 
            by a saturated or triply-unsaturated hydrocarbon chain (three isoprenoid units).

This implementation identifies the following features:
  1. A candidate 6-membered “dihydro” (pyran) ring with exactly one oxygen and five carbons that is not fully aromatic 
     (i.e. at least two atoms in the ring are non-aromatic).
  2. A candidate aromatic ring (benzene) defined as a 6‑membered ring with all atoms aromatic.
  3. A fused bicyclic system when the two rings share at least two atoms.
  4. A phenolic –OH substituent attached to the aromatic ring.
  5. An aliphatic side chain (only carbons, length ≥ 10 carbons) attached to an atom in the saturated ring that does not
     belong to the shared aromatic system. The side chain is required to be either fully saturated (0 double bonds) or contain 
     exactly 3 double bonds.
  
Note: This heuristic approach relaxes the previous strict exclusion of any aromatic atoms in the “saturated” ring.
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
    
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found in the molecule"
    
    # Find candidate saturated (dihydro) rings.
    # Look for 6-membered rings that contain exactly one oxygen and five carbons,
    # and not all the atoms are aromatic (i.e. at least 2 atoms are non-aromatic).
    candidate_saturated_rings = []
    for ring in ring_info:
        if len(ring) != 6:
            continue
        atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        atomic_nums = [atom.GetAtomicNum() for atom in atoms]
        if atomic_nums.count(8) != 1 or atomic_nums.count(6) != 5:
            continue
        # Allow rings that are not 100% aromatic.
        non_aromatic = sum(1 for atom in atoms if not atom.GetIsAromatic())
        if non_aromatic < 2:
            continue
        candidate_saturated_rings.append(set(ring))
    
    if not candidate_saturated_rings:
        return False, "No suitable 6-membered ring with one oxygen (dihydro part) found"
    
    # Find candidate aromatic rings (benzene-like).
    candidate_aromatic_rings = []
    for ring in ring_info:
        if len(ring) != 6:
            continue
        atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        if all(atom.GetIsAromatic() for atom in atoms):
            candidate_aromatic_rings.append(set(ring))
    
    if not candidate_aromatic_rings:
        return False, "No aromatic 6-membered ring (benzene moiety) found"
    
    # Identify a fused bicyclic core: one candidate saturated ring fused to an aromatic ring.
    found_core = False
    core_atoms = None
    core_saturated = None
    core_aromatic = None
    for sat_ring in candidate_saturated_rings:
        for aro_ring in candidate_aromatic_rings:
            if len(sat_ring.intersection(aro_ring)) >= 2:
                core_atoms = sat_ring.union(aro_ring)
                core_saturated = sat_ring
                core_aromatic = aro_ring
                found_core = True
                break
        if found_core:
            break
    if not found_core:
        return False, "No fused bicyclic chromanol core (dihydro ring fused to benzene) found"
    
    # Ensure the aromatic ring carries a phenolic -OH substituent.
    phenolic_found = False
    for idx in core_aromatic:
        atom = mol.GetAtomWithIdx(idx)
        # Check neighbors for an oxygen atom with at least one attached hydrogen.
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() > 0:
                phenolic_found = True
                break
        if phenolic_found:
            break
    if not phenolic_found:
        return False, "No phenolic hydroxyl (-OH) substituent found on the aromatic part of the core"
    
    # Identify the side chain attached at the 2‑position.
    # Strategy: in the saturated ring, use atoms that are not shared with the aromatic ring.
    side_chain_root = None
    sat_only = core_saturated.difference(core_aromatic)
    for idx in sat_only:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in core_atoms:
                continue
            if nbr.GetAtomicNum() == 6:  # expecting attachment via a carbon
                side_chain_root = nbr
                break
        if side_chain_root is not None:
            break
    if side_chain_root is None:
        return False, "No side chain attached to the non-aromatic portion of the saturated ring (expected at 2‑position)"
    
    # Extract the side chain using a breadth-first search from the attachment point,
    # excluding any atoms that belong to the core.
    side_chain_idxs = set()
    to_visit = [side_chain_root.GetIdx()]
    while to_visit:
        current = to_visit.pop()
        if current in side_chain_idxs:
            continue
        side_chain_idxs.add(current)
        for nbr in mol.GetAtomWithIdx(current).GetNeighbors():
            if nbr.GetIdx() in core_atoms:
                continue
            if nbr.GetIdx() not in side_chain_idxs:
                to_visit.append(nbr.GetIdx())
    
    if not side_chain_idxs:
        return False, "Failed to extract the side chain at the 2‑position"
    
    # Check that the side chain contains only carbon atoms.
    for idx in side_chain_idxs:
        if mol.GetAtomWithIdx(idx).GetAtomicNum() != 6:
            return False, f"Side chain contains a non-carbon atom: {mol.GetAtomWithIdx(idx).GetSymbol()} at index {idx}"
    
    n_chain_carbons = len(side_chain_idxs)
    if n_chain_carbons < 10:
        return False, f"Side chain is too short (found {n_chain_carbons} carbons; expected at least 10)"
    
    # Count double bonds present within the side chain.
    double_bond_count = 0
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        if a1 in side_chain_idxs and a2 in side_chain_idxs:
            if bond.GetBondType() == BondType.DOUBLE:
                double_bond_count += 1
    if double_bond_count not in (0, 3):
        return False, f"Side chain has {double_bond_count} double bond(s); expected 0 (saturated) or 3 (triply unsaturated)"
    
    return True, "Molecule contains a chromanol core (fused dihydropyran and benzene with phenolic -OH) with an appropriate 2‑position side chain"

# Testing examples (you may adjust or remove these tests as needed)
if __name__ == "__main__":
    test_smiles_list = [
        "CC(C)CCC[C@H](C)CCC[C@@H](C)CCC[C@@]1(C)CCC2=C(C)C(O)=C(C)C(C)=C2O1",  # (S,R,S)-alpha-tocopherol
        "CC1=C(C(=C2CCC(OC2=C1C)(C)CCCC(C)CCCC(C)CCCC(C)C)C)O",  # 2,5,7,8-tetramethyl-... (tocochromanol)
        "CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC[C@]1(C)CCc2c(C)c(O)c(C)c(C)c2O1",  # alpha-tocotrienol
        "CC(C)=CCC\\C(=C/CC\\C(C)=C\\CC[C@@]1(C)Oc2c(C)cc(O)cc2C=C1)C(O)=O",  # (R)-Sargachromenol
        "[H][C@@]1(CC\\C=C(/C)CC\\C=C(/C)CCC=C(C)C)CCc2cc(O)ccc2O1",  # didesmethyl tocotrienol
    ]
    for s in test_smiles_list:
        flag, reason = is_tocol(s)
        print("SMILES:", s)
        print("Result:", flag, "|", reason)
        print("------------------------------------------------")