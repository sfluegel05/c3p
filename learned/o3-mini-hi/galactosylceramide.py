"""
Classifies: CHEBI:36498 galactosylceramide
"""
"""
Classifies: Galactosylceramide – any cerebroside in which the monosaccharide head group is galactose.

A galactosylceramide is a glycosphingolipid consisting of a ceramide portion (a sphingoid base plus a fatty acid
attached via an amide bond) linked by an O‐glycosidic bond to a galactose sugar. 

This improved version applies several filters:
  1. The molecule is in the expected size range (30–150 heavy atoms).
  2. It has exactly one amide bond (the C(=O)N motif typical of ceramides).
  3. It contains at least one long, linear alkyl chain (here defined as a pattern of 8 consecutive carbons).
  4. It contains exactly one candidate sugar ring. Candidate rings are defined as six‐membered rings 
     that have exactly one ring oxygen and at least three exocyclic oxygen substituents.
     
This extra check on the number of candidate sugar rings is intended to filter out multi‐sugar glycosphingolipids
(such as glucosylceramides or disaccharide derivatives) that could otherwise be falsely flagged.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_galactosylceramide(smiles: str):
    """
    Determines if a molecule is a galactosylceramide based on its SMILES string.

    The algorithm applies several filters:
      1. Checks that the molecule is of moderate size (30–150 heavy atoms).
      2. Requires exactly one amide bond (the ceramide C(=O)N motif).
      3. Ensures the presence of at least one long alkyl chain (8 or more consecutive carbons).
      4. Searches for candidate sugar rings defined as six-membered rings with exactly one ring oxygen
         and at least three exocyclic oxygen substituents. Exactly one such candidate must be found.

    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a galactosylceramide, else False.
        str: Explanation for the classification decision.
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Check heavy atom count to ensure molecule is in a typical size range.
    heavy_atoms = mol.GetNumHeavyAtoms()
    if heavy_atoms < 30:
        return False, "Molecule appears too small to be a galactosylceramide."
    if heavy_atoms > 150:
        return False, "Molecule appears too large to be a typical galactosylceramide."
    
    # --- Step 1. Check for exactly one amide bond (the ceramide motif) ---
    # Using a SMARTS to capture a carbonyl (C=O) directly bonded to an N.
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) != 1:
        return False, f"Expected exactly 1 amide bond, found {len(amide_matches)}. (Not a typical ceramide)"
    
    # --- Step 2. Check for a long alkyl chain ---
    # A simple SMARTS pattern: 8 consecutive aliphatic carbons.
    chain_pattern = Chem.MolFromSmarts("CCCCCCCC")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "Missing a long alkyl chain (expected from fatty acyl or sphingoid base)."
    
    # --- Step 3. Look for candidate sugar rings ---
    # A candidate sugar ring is defined as one that:
    #    - has exactly 6 atoms (pyranose),
    #    - contains exactly one ring oxygen,
    #    - and has at least 3 exocyclic oxygen substituents on ring atoms.
    ring_info = mol.GetRingInfo()
    candidate_sugar_count = 0
    for ring in ring_info.AtomRings():
        if len(ring) != 6:
            continue  # only consider six-membered rings
        # Count ring oxygens.
        ring_oxygens = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        if ring_oxygens != 1:
            continue
        
        # Count external oxygen substituents (neighbors that are not in the ring).
        external_oxygen_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() == 8:
                    external_oxygen_count += 1
        if external_oxygen_count >= 3:
            candidate_sugar_count += 1

    if candidate_sugar_count == 0:
        return False, "No candidate pyranose sugar ring with sufficient oxygen substituents detected."
    if candidate_sugar_count > 1:
        return False, f"Found {candidate_sugar_count} candidate sugar rings; expected exactly 1 for a galactosylceramide."
    
    # --- All tests pass. Return a positive classification.
    return True, "Contains one ceramide amide bond, exactly one candidate pyranose sugar ring and a long alkyl chain, consistent with a galactosylceramide."

# Example usage:
# (Uncomment the following lines to test with a sample SMILES string.)
# smiles_example = "C(=C/CCCCCCCCCCCCC)\\[C@@H](O)[C@@H](NC(CCCCCCCCCCCCCCCCC)=O)CO[C@H]1[C@@H]([C@H]([C@H]([C@H](O1)CO)O)OS(O)(=O)=O)O"
# result, reason = is_galactosylceramide(smiles_example)
# print(result, reason)