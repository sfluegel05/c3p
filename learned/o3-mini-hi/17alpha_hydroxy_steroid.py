"""
Classifies: CHEBI:35342 17alpha-hydroxy steroid
"""
"""
Classifies: 17α-hydroxy steroid
Definition: The α-stereoisomer of 17-hydroxy steroid.

This version attempts to improve over the previous heuristic by:
  • Considering rings of size 5 or 6 as candidate parts of the steroid nucleus.
  • Clustering rings as “fused” if they share at least one atom (instead of 2) so as not to miss nuclei
    that are defined in nontraditional ways.
  • Accepting clusters where the union of atoms includes roughly 15–20 carbon atoms.
  • Requiring that the cluster contains at least one five‐membered ring (the presumptive D-ring).
  • Checking that at least one hydroxyl group (-OH) is attached to a carbon in the five‐membered ring
    that has a computed chiral descriptor (via the _CIPCode property).
    
Note: This heuristic can yield false positives or negatives because steroids are a very diverse class,
and our approach is a simplified approximation of the fused ring structural recognition.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_17alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is heuristically classified as a 17α-hydroxy steroid based on its SMILES.
    
    The method works by:
      1. Parsing the molecule and extracting ring information.
      2. Filtering to rings that are size 5 or 6 (the sizes typical in steroids).
      3. Clustering these rings into connected groups (if they share at least one atom).
      4. Looking for a cluster that has at least 2 (or, here, we require at least 3 for more robust steroid cores)
         rings and whose union contains between 15 and 20 carbon atoms.
      5. Requiring that one of the rings in that cluster is five-membered.
      6. Finally, checking that an -OH group (identified via a SMARTS query for O-hydrogen)
         is attached to a carbon in that five-membered ring that has a computed CIP stereochemistry label.
    
    Args:
       smiles (str): SMILES string of the molecule.
    
    Returns:
       bool: True if the molecule is classified as a 17α-hydroxy steroid, False otherwise.
       str: Explanation of the decision.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Get ring information from the molecule
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found in the molecule."
    
    # Filter to candidate rings of size 5 or 6 (only these sizes are used for the steroid core)
    candidate_rings = [set(ring) for ring in ring_info if len(ring) in (5, 6)]
    if not candidate_rings:
        return False, "No five- or six-membered rings found (cannot locate steroid nucleus)."
    
    # Cluster rings into connected groups using a more relaxed criteria: rings share at least 1 atom.
    clusters = []  # each cluster is a dict with 'rings' (list of sets) and 'union' (set of atom indices in cluster)
    for ring in candidate_rings:
        merged = False
        for cluster in clusters:
            # If the ring shares at least one atom with any ring in an existing cluster, merge it.
            if any(len(ring & existing) >= 1 for existing in cluster['rings']):
                cluster['rings'].append(ring)
                cluster['union'].update(ring)
                merged = True
                break
        if not merged:
            clusters.append({'rings': [ring], 'union': set(ring)})
    
    if not clusters:
        return False, "No fused ring clusters found."
    
    # Look among these clusters for one that can represent the steroid nucleus.
    steroid_cluster = None
    candidate_five_ring = None  # will store a five-membered ring within the cluster
    for cluster in clusters:
        # We require the cluster to have at least 2 rings; many steroids have four.
        if len(cluster['rings']) < 2:
            continue
        # Count carbon atoms in the cluster union:
        core_carbons = [idx for idx in cluster['union'] if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6]
        n_carbons = len(core_carbons)
        # Accept if number of carbon atoms is roughly in the range 15–20
        if 15 <= n_carbons <= 20:
            # Look for a five-membered ring in the cluster (a candidate D-ring)
            five_rings = [ring for ring in cluster['rings'] if len(ring) == 5]
            if five_rings:
                steroid_cluster = cluster
                candidate_five_ring = five_rings[0]  # choose the first five-membered ring as candidate
                break

    if steroid_cluster is None:
        msg = ("No fused ring cluster found with 2+ rings and 15–20 carbons in the union. "
               "The molecule deviates from the classical steroid nucleus.")
        return False, msg
    
    # Look for hydroxyl groups (-OH) using a SMARTS query for an O bonded to H.
    oh_query = Chem.MolFromSmarts("[OX2H]")
    oh_matches = mol.GetSubstructMatches(oh_query)
    if not oh_matches:
        return False, "No hydroxyl (-OH) group found in the molecule."
    
    # Now, check if at least one -OH is bonded to a carbon atom that is part of our candidate five-membered ring
    # and that the carbon has a computed CIP code (_CIPCode property) as a proxy for defined stereochemistry.
    hydroxyl_on_chiral = False
    for match in oh_matches:
        o_idx = match[0]
        o_atom = mol.GetAtomWithIdx(o_idx)
        # For each neighbor of the oxygen, check if it is a carbon in the candidate five-membered ring
        for nbr in o_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() in candidate_five_ring:
                # Check for computed CIP stereochemistry.
                if nbr.HasProp('_CIPCode'):
                    hydroxyl_on_chiral = True
                    break
        if hydroxyl_on_chiral:
            break

    if not hydroxyl_on_chiral:
        return False, ("No hydroxyl (-OH) group was found on a chiral carbon within the candidate D-ring "
                       "of the fused steroid core; cannot confirm 17α-hydroxy placement.")

    return True, ("Molecule's fused steroid core (cluster of %d rings with %d carbons overall, including at least one 5-membered ring) "
                  "and the presence of a hydroxyl on a chiral carbon within the five-membered ring heuristically classify it "
                  "as a 17α-hydroxy steroid." % (len(steroid_cluster['rings']), n_carbons))

# Example usage:
# Uncomment the following lines to test with one of the valid steroid SMILES:
# test_smiles = "[H][C@@]12CCC3=CC(=O)C=C[C@]3(C)[C@@]1([H])[C@@H](O)C[C@@]1(C)[C@@]2([H])C[C@@H](C(=O)OC)[C@]1(O)C(=O)CO"  # e.g., methyl prednisolone-16alpha-carboxylate
# result, reason = is_17alpha_hydroxy_steroid(test_smiles)
# print(result, reason)