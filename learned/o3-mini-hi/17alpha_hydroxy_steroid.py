"""
Classifies: CHEBI:35342 17alpha-hydroxy steroid
"""
"""
Classifies: 17α-hydroxy steroid
Definition: The α-stereoisomer of 17-hydroxy steroid.
This version attempts to improve on the previous heuristic by:
  • Considering rings of size 5 or 6 that are likely to belong to the steroid nucleus,
  • Clustering rings as fused if they share at least 2 atoms (rather than 1),
  • Looking for any fused cluster (of at least 2 rings) whose union has roughly 17 carbon atoms (16 to 18),
  • Requiring that one five‐membered ring is present in that cluster, and
  • Checking that at least one hydroxyl group (-OH) is attached to a chiral carbon (one with a CIP tag)
    within the five‐membered ring.
If these tests are passed, the molecule is heuristically classified as a 17α-hydroxy steroid.
Note: Because the stove‐pipe detection relies on multiple heuristics the outcome may occasionally be missed or be a false positive.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_17alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17α-hydroxy steroid based on its SMILES string.
    
    The method works by:
      - Parsing the molecule and extracting rings of size 5 or 6.
      - Clustering these rings as fused if they share at least 2 atoms.
      - Among clusters with at least 2 rings, finding one whose union (of atoms) contains 16–18 carbons.
      - Requiring that the cluster includes at least one five-membered ring.
      - Finally, checking that an -OH group (oxygen with a hydrogen, via SMARTS)
        is bonded to a carbon from the five-membered ring that also carries a CIP stereochemistry label.
        
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is heuristically classified as a 17α-hydroxy steroid,
              False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Get ring information (list of tuples of atom indices)
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found in the molecule."
    
    # Consider only rings that are size 5 or 6 (classical steroid rings).
    candidate_rings = [set(r) for r in ring_info if len(r) in (5, 6)]
    if not candidate_rings:
        return False, "No five- or six-membered rings found (cannot find steroid nucleus)."
    
    # Cluster rings into fused groups.
    # Here we consider two rings as fused if they share at least 2 atoms.
    clusters = []  # each cluster is a dict with keys: 'rings' (list of sets) and 'union' (set of atom indices)
    for ring in candidate_rings:
        merged = False
        for cluster in clusters:
            # Check if ring shares at least 2 atoms with any ring in this cluster.
            if any(len(ring & existing) >= 2 for existing in cluster['rings']):
                cluster['rings'].append(ring)
                cluster['union'].update(ring)
                merged = True
                break
        if not merged:
            clusters.append({'rings': [ring], 'union': set(ring)})
    
    if not clusters:
        return False, "No fused ring clusters found."
    
    # Now pick a candidate cluster that might be the steroid core.
    # We require at least two rings and that the union of atoms contains roughly 17 carbons.
    steroid_cluster = None
    for cluster in clusters:
        if len(cluster['rings']) < 2:
            continue
        # Count carbon atoms in the union
        core_carbons = [idx for idx in cluster['union'] if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6]
        n_carbons = len(core_carbons)
        if 16 <= n_carbons <= 18:
            # Also require that at least one ring in the cluster is five-membered.
            five_membered = [r for r in cluster['rings'] if len(r) == 5]
            if five_membered:
                steroid_cluster = cluster
                five_ring_atoms = five_membered[0]  # take the first 5-membered ring as candidate D-ring
                break

    if steroid_cluster is None:
        msg = "No fused ring cluster found with 2+ rings and ~17 carbons in the union. "
        msg += "This deviates from the classical steroid nucleus."
        return False, msg

    # Look for hydroxyl groups (-OH) using a SMARTS query.
    oh_query = Chem.MolFromSmarts("[OX2H]")
    oh_matches = mol.GetSubstructMatches(oh_query)
    if not oh_matches:
        return False, "No hydroxyl (-OH) group found in the molecule."
    
    # Check if at least one –OH is attached to a carbon in the candidate five-membered ring that has a CIP label.
    candidate_found = False
    for match in oh_matches:
        o_idx = match[0]
        o_atom = mol.GetAtomWithIdx(o_idx)
        # Look at all neighboring atoms
        for nbr in o_atom.GetNeighbors():
            # We want the neighbor to be a carbon that belongs to the five-membered ring
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() in five_ring_atoms:
                # Check if a CIP stereochemistry label has been computed.
                # (The presence of a '_CIPCode' property is a proxy for a defined stereochemistry.)
                if nbr.HasProp('_CIPCode'):
                    candidate_found = True
                    break
        if candidate_found:
            break
            
    if not candidate_found:
        return False, ("No hydroxyl (-OH) group found on a chiral carbon in the candidate five-membered (D-) ring; "
                       "cannot confirm 17α-hydroxy placement.")
    
    return True, ("Molecule's fused steroid core (a cluster of at least 2 fused rings with ~17 carbons, including one five-membered D-ring) "
                  "and the presence of a hydroxyl on a chiral carbon in the D-ring heuristically classify it as a 17α-hydroxy steroid.")

# Example usage (uncomment to test):
# test_smiles = "[H][C@@]12CCC3=CC(=O)C=C[C@]3(C)[C@@]1([H])[C@@H](O)C[C@@]1(C)[C@@]2([H])C[C@@H](C(=O)OC)[C@]1(O)C(=O)CO"  # methyl prednisolone-16alpha-carboxylate
# result, reason = is_17alpha_hydroxy_steroid(test_smiles)
# print(result, reason)