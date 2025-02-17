"""
Classifies: CHEBI:39437 tocol
"""
"""
Classifies: tocol
Definition: A chromanol with a chroman-6-ol skeleton that is substituted at position 2 
            by a saturated or triply-unsaturated (i.e. exactly 3 double bonds) hydrocarbon side chain 
            (consisting of three isoprenoid units).
            
Heuristic implementation:
  1. Identify a candidate dihydropyran (saturated) ring: a 6-membered ring with exactly one oxygen, 
     five carbons, and not fully aromatic.
  2. Identify a candidate aromatic (benzene-like) ring: a 6-membered ring with every atom aromatic.
  3. Find a fused bicyclic system if the saturated ring and aromatic ring share at least 2 atoms.
  4. Check that the aromatic ring carries a free (phenolic) –OH substituent.
  5. Identify the side chain attached at the “2‑position” on the dihydropyran using common chromanol connectivity:
     • In a proper chromanol the oxygen (in the saturated ring) is bound to two carbons. One of these is fused with 
       the aromatic ring while the other (non-fused) should carry the side chain.
  6. Extract the side chain (by breadth-first search from the attachment, excluding the core).
  7. Verify that the side chain is composed only of carbons, is long enough (≥10 carbons) and that its bonds 
     are either all single or (if unsaturated) contain exactly 3 double bonds.
     
Note: This heuristic does not fully resolve all edge cases but we have tried to improve the capture of the 2‑position 
side chain and the free –OH.
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
    
    # Get ring information
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found in the molecule"
    
    # Identify candidate saturated (dihydropyran) rings:
    cand_sat_rings = []
    for ring in ring_info:
        if len(ring) != 6:
            continue
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        # Should have exactly one oxygen and five carbons
        if sum(1 for a in atoms if a.GetAtomicNum() == 8) != 1 or sum(1 for a in atoms if a.GetAtomicNum() == 6) != 5:
            continue
        # Require that not all atoms are aromatic (i.e. at least 2 not aromatic)
        if sum(1 for a in atoms if not a.GetIsAromatic()) < 2:
            continue
        cand_sat_rings.append(set(ring))
    
    if not cand_sat_rings:
        return False, "No candidate dihydropyran (6-membered ring with one oxygen) found"
    
    # Identify candidate aromatic rings (benzene-like: 6-membered, all aromatic)
    cand_arom_rings = []
    for ring in ring_info:
        if len(ring) != 6:
            continue
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if all(a.GetIsAromatic() for a in atoms):
            cand_arom_rings.append(set(ring))
    
    if not cand_arom_rings:
        return False, "No candidate aromatic (benzene) ring found"
    
    # Look for a fused bicyclic system: any candidate saturated ring that shares at least 2 atoms with
    # any candidate aromatic ring.
    core_sat_ring = None
    core_arom_ring = None
    core_atoms = None
    for sat_ring in cand_sat_rings:
        for arom_ring in cand_arom_rings:
            if len(sat_ring.intersection(arom_ring)) >= 2:
                core_sat_ring = sat_ring
                core_arom_ring = arom_ring
                core_atoms = sat_ring.union(arom_ring)
                break
        if core_sat_ring is not None:
            break
    if core_sat_ring is None or core_arom_ring is None:
        return False, "No fused bicyclic chromanol core (dihydropyran fused to benzene) found"
    
    # Check that the aromatic ring carries a free (phenolic) -OH (an oxygen attached with at least one hydrogen)
    phenolic_found = False
    for idx in core_arom_ring:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() > 0:
                phenolic_found = True
                break
        if phenolic_found:
            break
    if not phenolic_found:
        return False, "No free phenolic -OH substituent found on the aromatic ring of the core"
    
    # Identify the side chain attachment at the expected 2‑position.
    # In the dihydropyran (core_sat_ring) identify the unique oxygen.
    oxygen_atom = None
    for idx in core_sat_ring:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 8:
            oxygen_atom = atom
            break
    if oxygen_atom is None:
        return False, "Core saturated ring does not contain an oxygen (unexpected)"
    
    # Get the neighbors of this oxygen that are in the saturated ring.
    pyran_neighbors = [nbr for nbr in oxygen_atom.GetNeighbors() if nbr.GetIdx() in core_sat_ring]
    if len(pyran_neighbors) != 2:
        return False, "Unexpected connectivity in dihydropyran: oxygen does not have two ring neighbors"
    
    # One neighbor should be part of the fused aromatic ring and the other not.
    pos2_atom = None
    for nbr in pyran_neighbors:
        if nbr.GetIdx() not in core_arom_ring:
            pos2_atom = nbr
            break
    if pos2_atom is None:
        return False, "Could not determine the non-fused neighbor of the oxygen (expected 2‑position)"
    
    # From the candidate “2‑position” atom, find a substituent that is not part of the core.
    side_chain_root = None
    for nbr in pos2_atom.GetNeighbors():
        if nbr.GetIdx() in core_atoms:
            continue
        if nbr.GetAtomicNum() == 6:  # should be a carbon
            side_chain_root = nbr
            break
    if side_chain_root is None:
        return False, "No side chain attached at the expected 2‑position of the chromanol core"
        
    # Extract the side chain: BFS from the side_chain_root excluding any atom in the core.
    side_chain_idxs = set()
    to_visit = [side_chain_root.GetIdx()]
    while to_visit:
        current = to_visit.pop()
        if current in side_chain_idxs:
            continue
        side_chain_idxs.add(current)
        atom = mol.GetAtomWithIdx(current)
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in core_atoms:
                continue
            if nbr.GetIdx() not in side_chain_idxs:
                to_visit.append(nbr.GetIdx())
                
    if not side_chain_idxs:
        return False, "Failed to extract the side chain at the 2‑position"
    
    # Check side chain is entirely carbon atoms.
    for idx in side_chain_idxs:
        if mol.GetAtomWithIdx(idx).GetAtomicNum() != 6:
            return False, f"Side chain contains a non-carbon atom: {mol.GetAtomWithIdx(idx).GetSymbol()} at index {idx}"
    
    # Verify chain length (number of carbons): a chain of three isoprenoid units is usually ≥10.
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
        return False, f"Side chain has {double_bond_count} double bond(s); should be either fully saturated (0) or with exactly 3 double bonds"
    
    return True, "Molecule contains a chromanol core (dihydropyran fused to benzene with free -OH) with an appropriate 2‑position side chain"

# Testing examples (you may adjust or remove tests as needed)
if __name__ == "__main__":
    test_smiles_list = [
        "CC(C)CCC[C@H](C)CCC[C@@H](C)CCC[C@@]1(C)CCC2=C(C)C(O)=C(C)C(C)=C2O1",  # (S,R,S)-alpha-tocopherol
        "CC1=C(C(=C2CCC(OC2=C1C)(C)CCCC(C)CCCC(C)CCCC(C)C)C)O",  # 2,5,7,8-tetramethyl-... (side chain may be short)
        "CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC[C@]1(C)CCc2c(C)c(O)c(C)c(C)c2O1",  # alpha-tocotrienol
        "CC(C)=CCC\\C(=C/CC\\C(C)=C\\CC[C@@]1(C)Oc2c(C)cc(O)cc2C=C1)C(O)=O",  # (R)-Sargachromenol (side chain contains O, invalid)
        "[H][C@@]1(CC\\C=C(/C)CC\\C=C(/C)CCC=C(C)C)CCc2cc(O)ccc2O1",  # didesmethyl tocotrienol
    ]
    for s in test_smiles_list:
        flag, reason = is_tocol(s)
        print("SMILES:", s)
        print("Result:", flag, "|", reason)
        print("------------------------------------------------")