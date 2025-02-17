"""
Classifies: CHEBI:1722 3beta-hydroxy-Delta(5)-steroid
"""
"""
Classifies: 3beta-hydroxy-Delta(5)-steroid
Definition: Any 3β-hydroxy-steroid that contains a double bond between positions 5 and 6.
Improvement: In addition to requiring a fused steroid nucleus (≥3 rings), we now look for a 
six-membered ring that contains both a double bond and, in a position that is not directly 
involved in the double bond, a chiral (stereogenic) carbon with an –OH substituent. We also 
require that the –OH group be at least two atoms away (along the ring) from the double bond.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3beta_hydroxy_Delta_5__steroid(smiles: str):
    """
    Determines if the given SMILES string represents a 3β-hydroxy-Δ5-steroid.

    This function applies several heuristics:
      1. It extracts all rings from the molecule and builds a graph of fused rings (rings that share ≥2 atoms).
         It then selects the largest connected (fused) ring system. A steroid nucleus is assumed to have ≥3 rings.
      2. It iterates over each six-membered ring in the nucleus. For each ring it checks:
           (a) Does the ring contain a double bond? (We require at least one bond in the ring to be double.)
           (b) Does the ring contain a chiral (stereogenic) carbon with at least one hydroxyl (–OH) substituent?
               Moreover, the candidate –OH bearing carbon should not be one of the atoms forming the double bond,
               and (as a crude proxy for position 3) its distance (in ring order) from either end of the double bond should be >= 2.
      3. Only if at least one six‐membered ring in the steroid nucleus meets these requirements do we return True.

    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets the criteria, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES into an RDKit Mol
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information (each ring is a tuple of atom indices)
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings detected"
    
    # Build an adjacency for rings: two rings are fused if they share at least 2 atoms.
    n_rings = len(ring_info)
    ring_adj = {i: set() for i in range(n_rings)}
    for i in range(n_rings):
        atoms_i = set(ring_info[i])
        for j in range(i+1, n_rings):
            atoms_j = set(ring_info[j])
            if len(atoms_i.intersection(atoms_j)) >= 2:
                ring_adj[i].add(j)
                ring_adj[j].add(i)
    
    # Find clusters (connected components) of fused rings using DFS.
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
    
    # Get a set of all atom indices belonging to the nucleus.
    nucleus_atoms = set()
    for ring_idx in largest_cluster:
        nucleus_atoms.update(ring_info[ring_idx])
    
    # Look for a six-membered ring in the nucleus that meets the Δ5 and 3β-OH criteria.
    for ring_idx in largest_cluster:
        ring = list(ring_info[ring_idx])
        if len(ring) != 6:
            continue  # must be a six-membered ring
        
        # Identify bonds within the ring that are double bonds.
        # Also record the indices (within the ring ordering) of the atoms involved.
        double_bond_positions = []
        ring_order = ring[:]  # the order as provided by RDKit ring info
        n = len(ring_order)
        for i in range(n):
            a1 = ring_order[i]
            a2 = ring_order[(i+1) % n]
            bond = mol.GetBondBetweenAtoms(a1, a2)
            if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                double_bond_positions.append(i)  # record position of the first atom of this bond
        
        if not double_bond_positions:
            continue  # no double bond in this ring
        
        # Now check for a chiral carbon with a hydroxyl group in the same ring.
        # The candidate should not be one of the double bond atoms and should be at least 2 bonds away (along the ring)
        # from any atom that is part of a double bond.
        for pos, atom_idx in enumerate(ring_order):
            atom = mol.GetAtomWithIdx(atom_idx)
            # Skip if this atom is involved in any double bond in the ring.
            involved_in_db = False
            for db_pos in double_bond_positions:
                # double bond spans positions db_pos and db_pos+1 mod n
                if pos == db_pos or pos == (db_pos+1)%n:
                    involved_in_db = True
                    break
            if involved_in_db:
                continue
            
            # Check if atom has a chiral tag (_CIPCode property can be used as a proxy for stereochemistry)
            if not atom.HasProp("_CIPCode"):
                continue

            # Check if atom is sp3 carbon
            if atom.GetAtomicNum() != 6 or atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                continue

            # Look among neighbors for an oxygen that is part of an –OH group.
            # We require that the oxygen bears at least one hydrogen.
            has_OH = False
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8:
                    # Check bonded hydrogens (include implicit and explicit)
                    if nbr.GetTotalNumHs() >= 1:
                        has_OH = True
                        break
            if not has_OH:
                continue

            # Check ring-distance from this candidate atom to each double bond position.
            valid_distance = False
            for db_pos in double_bond_positions:
                # Calculate circular distance along the ring.
                d1 = abs(pos - db_pos)
                d2 = n - d1
                if min(d1, d2) >= 2:
                    valid_distance = True
                    break
            if not valid_distance:
                continue
            
            # Found a valid candidate in this ring!
            return True, ("Molecule contains a fused steroid nucleus (≥3 fused rings) with a six-membered ring "
                          "that has a double bond and a chiral carbon bearing an –OH group (not part of the double bond), "
                          "consistent with a 3β-hydroxy-Δ5-steroid.")
    
    return False, "No six-membered ring in the steroid nucleus with both a double bond and an appropriately positioned hydroxyl group was found"