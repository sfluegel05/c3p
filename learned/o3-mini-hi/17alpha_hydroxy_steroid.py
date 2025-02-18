"""
Classifies: CHEBI:35342 17alpha-hydroxy steroid
"""
"""
Classifies: 17α-hydroxy steroid
Definition: The α-stereoisomer of 17-hydroxy steroid.
This code attempts to detect a classical steroid nucleus by:
  • Considering only rings of size 5 or 6 in the full molecule,
  • Clustering rings that share at least one atom (“fused” rings),
  • Requiring that the best fused cluster has at least 4 rings,
  • Checking that the union of atoms (counting only carbons) in that cluster is exactly 17,
  • Requiring that exactly one five‐membered ring is present (the D‐ring), and
  • Verifying that at least one hydroxyl group (–OH) lies on a chiral carbon within a five‐membered ring.
If these tests are passed, the molecule is heuristically classified as a 17α-hydroxy steroid.
Note: This approach is heuristic and may not capture every nuance.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_17alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17α-hydroxy steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a 17α-hydroxy steroid, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Get all ring information from the original molecule.
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found in the molecule."

    # Filter to rings of size 5 or 6 (typical of classical steroid systems).
    candidate_rings = [set(r) for r in ring_info if len(r) in (5, 6)]
    if not candidate_rings:
        return False, "No five- or six-membered rings found (cannot find steroid nucleus)."
        
    # Cluster rings that are “fused”. Here we use a simple criterion:
    # two rings are considered fused if they share at least one atom.
    clusters = []  # each cluster is a dict with keys 'rings' (list of sets) and 'union' (set of atom indices)
    for ring in candidate_rings:
        merged = False
        for cluster in clusters:
            # if the new ring shares at least one atom with any ring already in the cluster, merge it.
            if any(len(ring & existing) >= 1 for existing in cluster['rings']):
                cluster['rings'].append(ring)
                cluster['union'].update(ring)
                merged = True
                break
        if not merged:
            clusters.append({'rings': [ring], 'union': set(ring)})
    
    if not clusters:
        return False, "No fused ring clusters found."
    
    # Choose the cluster with the greatest number of rings (heuristically, the steroid core)
    best_cluster = max(clusters, key=lambda c: len(c['rings']))
    
    fused_ring_count = len(best_cluster['rings'])
    if fused_ring_count < 4:
        return False, f"Stereocore not detected: Found only {fused_ring_count} fused rings; expected at least 4 for a classical steroid nucleus."
    
    # Count the number of carbon atoms in the union of all fused rings.
    core_atom_indices = best_cluster['union']
    core_carbons = [idx for idx in core_atom_indices if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6]
    if len(core_carbons) != 17:
        return False, f"Fused core does not have 17 carbons (found {len(core_carbons)}). This deviates from the classical steroid nucleus."
    
    # Check for the presence of a five-membered ring in the cluster.
    # Classical steroids have one D-ring of size 5.
    five_membered_rings = [r for r in best_cluster['rings'] if len(r) == 5]
    if len(five_membered_rings) != 1:
        return False, f"Expected exactly one five-membered ring in the steroid core, found {len(five_membered_rings)}."
    
    # Identify the atoms in the single five-membered ring.
    five_ring_atoms = five_membered_rings[0]
    
    # Look for hydroxyl groups (-OH) in the molecule.
    oh_query = Chem.MolFromSmarts("[OX2H]")
    oh_matches = mol.GetSubstructMatches(oh_query)
    if not oh_matches:
        return False, "No hydroxyl (-OH) group found in the molecule."

    # Check if at least one –OH is attached to a carbon (atomic number 6)
    # that is part of the five-membered ring and which has an assigned CIP stereochemistry.
    candidate_found = False
    for match in oh_matches:
        o_idx = match[0]  # the oxygen atom index in the OH group
        o_atom = mol.GetAtomWithIdx(o_idx)
        for nbr in o_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() in five_ring_atoms:
                # Check that a CIP label has been assigned to this carbon.
                if nbr.HasProp('_CIPCode'):
                    candidate_found = True
                    break
        if candidate_found:
            break
    if not candidate_found:
        return False, ("No hydroxyl group found on a chiral carbon in the fused five-membered (D-) ring; "
                       "cannot confirm 17α-hydroxy placement.")
    
    return True, ("Molecule's fused steroid core (four fused rings with 17 carbons including a five-membered D-ring) "
                  "and the presence of a hydroxyl on a chiral carbon in the D-ring heuristically classify it as a 17α-hydroxy steroid.")

# Example usage:
# test_smiles = "[H][C@@]12CCC3=CC(=O)C=C[C@]3(C)[C@@]1([H])[C@@H](O)C[C@@]1(C)[C@@]2([H])C[C@@H](C(=O)OC)[C@]1(O)C(=O)CO"  # methyl prednisolone-16alpha-carboxylate
# result, reason = is_17alpha_hydroxy_steroid(test_smiles)
# print(result, reason)