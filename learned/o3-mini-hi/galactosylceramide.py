"""
Classifies: CHEBI:36498 galactosylceramide
"""
"""
Classifies: Galactosylceramide – any cerebroside in which the monosaccharide head group is galactose.

A galactosylceramide is a glycosphingolipid consisting of a ceramide (sphingoid base plus fatty acid, attached via an amide bond)
linked by an O‐glycosidic bond to a galactose sugar. In this version the molecule must satisfy:
  1. Size between 30 and 150 heavy atoms.
  2. Exactly one amide bond (C(=O)N) typical for the ceramide.
  3. At least one long (8+ carbon) uninterrupted alkyl chain.
  4. Exactly one candidate 6‐membered “sugar ring” (a pyranose) that has one ring oxygen and at least three exocyclic oxygens.
  5. The candidate sugar ring must further match a galactose motif (it must have the stereochemical “fingerprint”
     of either α‐ or β‐galactopyranose).
     
Because the latter step directly checks that the headgroup is galactose (and not, say, glucose) the false positives are reduced.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_galactosylceramide(smiles: str):
    """
    Determines if a molecule is a galactosylceramide based on its SMILES string.

    The algorithm applies several filters:
      1. Molecule size is within a typical range for glycosphingolipids.
      2. Exactly one amide bond (the ceramide motif).
      3. Presence of at least one long aliphatic chain (≥8 continuous carbons).
      4. Exactly one candidate six‐membered sugar ring (with one ring O and at least 3 exocyclic O substituents).
      5. The candidate sugar ring must match a galactose pattern (either alpha or beta).

    Args:
      smiles (str): SMILES string of the molecule.
        
    Returns:
      bool: True if the molecule is classified as a galactosylceramide.
      str: Explanation for the decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # --- Filter 1. Check heavy atom count (typical sugars + ceramide lie in a moderate size range) ---
    heavy_atoms = mol.GetNumHeavyAtoms()
    if heavy_atoms < 30:
        return False, "Molecule appears too small to be a galactosylceramide."
    if heavy_atoms > 150:
        return False, "Molecule appears too large to be a typical galactosylceramide."
    
    # --- Filter 2. Check for exactly one amide bond (the ceramide motif) ---
    # SMARTS: a carbonyl directly bonded to a nitrogen.
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) != 1:
        return False, f"Expected exactly 1 amide bond, found {len(amide_matches)} (not a typical ceramide)."
    
    # --- Filter 3. Check for a long aliphatic chain (≥8 consecutive carbons) ---
    chain_pattern = Chem.MolFromSmarts("CCCCCCCC")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "Missing a long alkyl chain (expected from a fatty acyl or sphingoid base)."
    
    # --- Filter 4. Find candidate sugar rings ---
    # A candidate sugar ring is defined by:
    #   - Exactly 6 atoms in the ring,
    #   - Exactly one of those atoms is oxygen,
    #   - At least 3 oxygen substituents on ring atoms (outside the ring).
    ring_info = mol.GetRingInfo()
    candidate_ring_indices = []
    for ring in ring_info.AtomRings():
        if len(ring) != 6:
            continue  # not a pyranose ring
        # Count number of ring oxygens.
        n_ring_oxygens = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        if n_ring_oxygens != 1:
            continue  # not the right type of ring
        # Count exocyclic oxygen substituents on atoms of the ring.
        ext_oxygen_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() == 8:
                    ext_oxygen_count += 1
        if ext_oxygen_count >= 3:
            candidate_ring_indices.append(ring)
    
    if len(candidate_ring_indices) == 0:
        return False, "No candidate pyranose sugar ring with sufficient oxygen substituents detected."
    if len(candidate_ring_indices) > 1:
        return False, f"Found {len(candidate_ring_indices)} candidate sugar rings; expected exactly 1 for a galactosylceramide."
    
    # --- Filter 5. Verify that the candidate sugar ring is galactose.
    # We create two SMARTS queries for the galactose motif.
    # Note: These patterns assume a pyranose ring (6-membered) with five hydroxyl groups.
    # They differ in the stereochemistry at the ring’s anomeric (linkage) carbon.
    alpha_gal_smarts = "OC[C@H]1O[C@@H](O)[C@@H](O)[C@H](O)[C@@H]1O"
    beta_gal_smarts  = "OC[C@@H]1O[C@@H](O)[C@@H](O)[C@H](O)[C@@H]1O"
    alpha_gal = Chem.MolFromSmarts(alpha_gal_smarts)
    beta_gal  = Chem.MolFromSmarts(beta_gal_smarts)
    
    # Extract the candidate sugar ring as a submol.
    # (Use the indices of the ring; note that the resulting submol does not include attachments.)
    candidate_ring = candidate_ring_indices[0]
    ring_mol = Chem.PathToSubmol(mol, candidate_ring)
    
    # Now test if the extracted ring (or parts of it) match one of our galactose SMARTS patterns.
    # We are willing to accept either an alpha or beta configuration.
    if not (ring_mol.HasSubstructMatch(alpha_gal) or ring_mol.HasSubstructMatch(beta_gal)):
        return False, "Candidate sugar ring does not display the expected galactose stereochemistry."
    
    # All tests pass.
    return True, "Contains one ceramide amide bond, exactly one candidate pyranose sugar ring with galactose stereochemistry, and a long alkyl chain – consistent with a galactosylceramide."

# Example usage:
# Uncomment the lines below to test one of the provided SMILES:
# smiles_example = "C(=C/CCCCCCCCCCCCC)\\[C@@H](O)[C@@H](NC(CCCCCCCCCCCCCCCCC)=O)CO[C@H]1[C@@H]([C@H]([C@H]([C@H](O1)CO)O)OS(O)(=O)=O)O"
# result, reason = is_galactosylceramide(smiles_example)
# print(result, reason)