"""
Classifies: CHEBI:1722 3beta-hydroxy-Delta(5)-steroid
"""
"""
Classifies: 3beta-hydroxy-Delta(5)-steroid
Definition: Any 3β-hydroxy-steroid that contains a double bond between positions 5 and 6.
Improvement: Force stereochemistry assignment and kekulization, then require 
(a) a fused ring nucleus (≥3 rings), and (b) at least one six-membered ring within that nucleus 
has a nonaromatic double bond and a chiral (sp³) carbon bearing a free –OH group; further, 
this –OH-bearing carbon must be at least two atoms away (circulatory) from the ends of any double bond.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3beta_hydroxy_Delta_5__steroid(smiles: str):
    """
    Checks whether the SMILES corresponds to a 3β-hydroxy-Δ5-steroid.
    
    The algorithm works in several steps:
    1. Parse the SMILES and assign stereochemistry; kekulize to make double bonds explicit.
    2. Get all rings and group fused ones (rings sharing ≥2 atoms). The largest fused ring system 
       is taken as a steroid nucleus candidate. (A steroid nucleus is assumed to have ≥3 rings.)
    3. For each six-membered ring in the nucleus, we first locate any nonaromatic double bonds.
       Then, we examine each atom in the ring. For an atom to be a candidate for the 3β–OH:
         (a) it must be an sp³ carbon and show chiral information (via its chiral tag).
         (b) It should have at least one neighbor oxygen that is “free” (i.e. bound to hydrogen).
         (c) It must not be part of any double bond, and its position in the ring must be at least
             2 steps away (circularly) from every atom that participates in a double bond.
    4. If at least one six‐membered ring meets these conditions, return True.
    
    Args:
        smiles (str): The SMILES string.
    
    Returns:
        bool: True if the molecule is classified as a 3β-hydroxy-Δ5-steroid, False otherwise.
        str: Explanation for the classification.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception:
        return False, "Error parsing SMILES"
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure stereochemistry is assigned and kekulize the molecule so double bonds are explicit.
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception:
        # Some molecules may not kekulize; we continue nonetheless.
        pass

    # Get ring information (each ring is a tuple of atom indices)
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings detected"
    
    # Build fused ring connectivity: two rings are fused if they share at least 2 atoms.
    n_rings = len(ring_info)
    ring_adj = {i: set() for i in range(n_rings)}
    for i in range(n_rings):
        set_i = set(ring_info[i])
        for j in range(i+1, n_rings):
            set_j = set(ring_info[j])
            if len(set_i.intersection(set_j)) >= 2:
                ring_adj[i].add(j)
                ring_adj[j].add(i)
    
    # Find connected clusters of fused rings via DFS.
    visited = set()
    clusters = []
    for i in range(n_rings):
        if i not in visited:
            cluster = set()
            stack = [i]
            while stack:
                cur = stack.pop()
                if cur in cluster:
                    continue
                cluster.add(cur)
                stack.extend(ring_adj[cur] - cluster)
            clusters.append(cluster)
            visited |= cluster

    # Choose the largest fused ring cluster as the steroid nucleus candidate.
    largest_cluster = max(clusters, key=lambda x: len(x))
    if len(largest_cluster) < 3:
        return False, "No fused ring system (steroid nucleus) of at least 3 rings found"
    
    # Get all atom indices that are part of this nucleus.
    nucleus_atoms = set()
    for ring_idx in largest_cluster:
        nucleus_atoms.update(ring_info[ring_idx])
    
    # Loop over each ring in the nucleus candidate that is a six-membered ring.
    for ring_idx in largest_cluster:
        ring = list(ring_info[ring_idx])
        if len(ring) != 6:
            continue  # Only six-membered rings are candidates
        
        # Identify bonds within the ring that are nonaromatic double bonds.
        double_bond_positions = []
        n = len(ring)
        # Use the order as provided by the ring (assumed cyclic)
        ring_order = ring[:]  
        for i in range(n):
            a1 = ring_order[i]
            a2 = ring_order[(i+1) % n]
            bond = mol.GetBondBetweenAtoms(a1, a2)
            if bond is None:
                continue
            # Only count if the bond is a double bond and not flagged aromatic.
            if bond.GetBondType() == Chem.BondType.DOUBLE and not bond.GetIsAromatic():
                double_bond_positions.append(i)  # record the index (first atom in bond)
        
        if not double_bond_positions:
            continue  # this ring does not have the needed double bond
        
        # Look for a candidate atom in the ring that carries an –OH group.
        for pos, atom_idx in enumerate(ring_order):
            atom = mol.GetAtomWithIdx(atom_idx)
            # Skip if atom is part of any double bond in the ring.
            in_double = False
            for db_pos in double_bond_positions:
                # double bond spans index db_pos and (db_pos+1)%n
                if pos == db_pos or pos == ((db_pos+1) % n):
                    in_double = True
                    break
            if in_double:
                continue
            
            # Check that the atom is carbon, sp3 and marked as chiral by its chiral tag.
            if atom.GetAtomicNum() != 6:
                continue
            if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                continue
            if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
                continue

            # Look for an –OH group attached to this atom.
            has_OH = False
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8:
                    # Require that oxygen has at least one hydrogen (free hydroxyl group).
                    if nbr.GetTotalNumHs() >= 1:
                        has_OH = True
                        break
            if not has_OH:
                continue
            
            # Check that this atom is at least 2 bonds (cyclically) away from every atom involved in a double bond.
            valid_distance = True
            for db_pos in double_bond_positions:
                # For each double bond, check both positions.
                for db_atom in [db_pos, (db_pos+1) % n]:
                    d1 = abs(pos - db_atom)
                    d2 = n - d1
                    if min(d1, d2) < 2:
                        valid_distance = False
                        break
                if not valid_distance:
                    break
            if not valid_distance:
                continue
            
            # If we reach here, we found a candidate six-membered ring that meets the criteria.
            return True, ("Molecule contains a fused steroid nucleus (≥3 fused rings) and a six-membered ring "
                          "with a nonaromatic double bond plus a chiral sp³ carbon bearing a free –OH group located "
                          "at least 2 bonds away (along the ring) from the double bond. This is consistent with a "
                          "3β-hydroxy-Δ5-steroid.")
    
    return False, ("No six-membered ring in the steroid nucleus was found that simultaneously has a "
                   "nonaromatic double bond and an appropriately positioned (chiral, OH-bearing) carbon.")
    
# Example usage:
if __name__ == "__main__":
    # Here one may test a few examples (this is not exhaustive)
    test_smiles = ["C[C@@H](CO)[C@@H](O)CC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C",  # (24S,25S)-cholest-5-en-3beta,24,26-triol
                   "CC[C@H](C)C(=O)[C@@H]1[C@H]2[C@@H]3[C@@H](O)OC(C)=CC3=CC(=O)[C@]2(C)OC1=O"]  # a false positive example
    for s in test_smiles:
        flag, reason = is_3beta_hydroxy_Delta_5__steroid(s)
        print("SMILES:", s)
        print("Classification:", flag)
        print("Reason:", reason)
        print("------------------------------------------------")