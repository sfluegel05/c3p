"""
Classifies: CHEBI:35342 17alpha-hydroxy steroid
"""
"""
Classifies: 17α-hydroxy steroid
Definition: The α-stereoisomer of 17-hydroxy steroid.
This improved code uses a fused–ring clustering approach on the Murcko scaffold
to detect a classical steroid core. It then checks that:
  • The fused core comprises at least four rings (with at least one five‐member ring),
  • The union of atoms in the core has exactly 17 carbon atoms,
  • At least one hydroxyl (–OH) group is attached to a chiral carbon in a five‐member ring.
If these tests are passed, the molecule is heuristically classified as a 17α-hydroxy steroid.
Note: A complete stereochemical assignment at C17 is not performed.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_17alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17α-hydroxy steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a 17α-hydroxy steroid, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Obtain the Murcko scaffold (the core framework) of the molecule.
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold is None:
        return False, "Could not obtain Murcko scaffold."

    # Get ring information from the scaffold.
    ring_info = scaffold.GetRingInfo()
    rings = [set(r) for r in ring_info.AtomRings()]
    if not rings:
        return False, "No rings found in the scaffold."

    # Cluster the rings that are fused (i.e. share at least one atom)
    clusters = []
    for r in rings:
        merged = False
        for cluster in clusters:
            if cluster & r:  # if there is an intersection between ring r and the cluster
                cluster.update(r)
                merged = True
                break
        if not merged:
            clusters.append(set(r))
    # For our purposes choose the largest fused set
    core_cluster = max(clusters, key=lambda c: len(c))
    
    # In addition, count how many individual rings (from the original rings list)
    # have all their atoms included in the core cluster.
    fused_ring_count = 0
    five_membered_present = False
    for r in rings:
        if r.issubset(core_cluster):
            fused_ring_count += 1
            if len(r) == 5:
                five_membered_present = True
    if fused_ring_count < 4:
        return False, (f"Stereocore not detected: Found only {fused_ring_count} fused rings in the scaffold; "
                       "expected at least 4 for a classical steroid nucleus.")
    if not five_membered_present:
        return False, "No five-membered ring (D-ring candidate) detected in the fused core."

    # Count the number of carbon atoms in the core cluster.
    # (Note: these indices are those of the scaffold molecule.)
    core_carbon_count = sum(1 for idx in core_cluster if scaffold.GetAtomWithIdx(idx).GetAtomicNum() == 6)
    if core_carbon_count != 17:
        return False, (f"Fused core does not have 17 carbons (found {core_carbon_count}). "
                       "This deviates from the classical steroid nucleus.")

    # Now, to verify the hydroxyl position we need to map the scaffold core back to the original molecule.
    scaffold_matches = mol.GetSubstructMatches(scaffold)
    if not scaffold_matches:
        return False, "Unable to map the scaffold back to the original molecule."
    # We use the first mapping. mapping[i] gives the atom index in the molecule corresponding 
    # to atom index i in the scaffold.
    mapping = scaffold_matches[0]
    
    # Identify the atom indices in the original molecule that correspond to atoms in the five-membered ring(s)
    five_ring_scaffold_indices = set()
    for r in rings:
        if len(r) == 5 and r.issubset(core_cluster):
            five_ring_scaffold_indices.update(r)
    five_ring_mol_indices = set(mapping[i] for i in five_ring_scaffold_indices)

    # Look for hydroxyl groups in the molecule using a simple SMARTS for –OH.
    oh_query = Chem.MolFromSmarts("[OX2H]")
    oh_matches = mol.GetSubstructMatches(oh_query)
    if not oh_matches:
        return False, "No hydroxyl (-OH) group found in the molecule."

    # Check if at least one –OH group is attached to a carbon (atomic number 6) that is:
    #   (a) part of the fused five-membered ring in the core, and 
    #   (b) is chiral (has a CIP label).
    candidate_found = False
    for match in oh_matches:
        o_idx = match[0]
        o_atom = mol.GetAtomWithIdx(o_idx)
        for neighbor in o_atom.GetNeighbors():
            # Look for a carbon that is in the five-membered core
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() in five_ring_mol_indices:
                # Check for assigned CIP stereochemistry on this carbon.
                if neighbor.HasProp('_CIPCode'):
                    candidate_found = True
                    break
        if candidate_found:
            break
    if not candidate_found:
        return False, ("No hydroxyl group found on a chiral carbon in the fused five-membered ring "
                       "of the core; cannot confirm 17α-hydroxy placement.")

    return True, ("Molecule's fused steroid core (detected as at least 4 fused rings with 17 carbons, including "
                  "a five-membered D-ring) and the presence of a hydroxyl on a chiral carbon in the D-ring "
                  "heuristically classify it as a 17α-hydroxy steroid.")

# Example usage:
# test_smiles = "[H][C@@]12CCC3=CC(=O)C=C[C@]3(C)[C@@]1([H])[C@@H](O)C[C@@]1(C)[C@@]2([H])C[C@@H](C(=O)OC)[C@]1(O)C(=O)CO"  # methyl prednisolone-16alpha-carboxylate
# result, reason = is_17alpha_hydroxy_steroid(test_smiles)
# print(result, reason)