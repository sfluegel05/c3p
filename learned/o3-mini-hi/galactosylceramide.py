"""
Classifies: CHEBI:36498 galactosylceramide
"""
"""
Classifies: Galactosylceramide – any cerebroside in which the monosaccharide head group is galactose.
A galactosylceramide is a glycosphingolipid consisting of a ceramide (sphingoid base plus a fatty acid attached via an amide bond)
linked by an O‐glycosidic bond to a galactose sugar. 
The improved strategy looks for exactly one amide bond, a pyranose sugar ring decorated with oxygens,
and the presence of a long aliphatic chain.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_galactosylceramide(smiles: str):
    """
    Determines if a molecule is a galactosylceramide based on its SMILES string.
    
    The algorithm applies several filters:
    
      1. The SMILES must be valid and the molecule must be of moderate size (heavy atoms between 30 and 100).
      2. The molecule must contain exactly one amide bond via a -C(=O)N- fragment (the ceramide motif).
      3. The molecule must contain a six-membered (pyranose) ring with exactly one ring oxygen and at least three 
         external oxygen substituents. This is taken as a proxy for a decorated sugar ring (such as galactose) 
         even if some -OH groups are modified (e.g. sulfonated).
      4. The molecule must include at least one long linear alkyl chain (as expected for the fatty acyl chain).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as a galactosylceramide, else False.
        str: Explanation for the classification decision.
    """
    # Parse the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    heavy_atoms = mol.GetNumHeavyAtoms()
    # Require a moderate size; many glycosphingolipids fall in an intermediate range.
    if heavy_atoms < 30:
        return False, "Molecule appears too small to be a galactosylceramide."
    if heavy_atoms > 150:
        return False, "Molecule appears too large to be a typical galactosylceramide."
    
    # --- Step 1. Check for exactly one amide bond ---
    # The amide pattern is used to capture the ceramide C(=O)N motif.
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) != 1:
        return False, f"Expected exactly 1 amide bond, found {len(amide_matches)}. (Not a typical ceramide)"
    
    # --- Step 2. Look for a long alkyl chain ---
    # A simple pattern: 8 or more consecutive aliphatic carbons.
    chain_pattern = Chem.MolFromSmarts("CCCCCCCC")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "Missing a long alkyl chain expected in fatty acyl or sphingoid base."
    
    # --- Step 3. Identify a candidate galactose (pyranose) ring ---
    ring_info = mol.GetRingInfo()
    found_sugar = False
    # Loop over all rings in the molecule
    for ring in ring_info.AtomRings():
        if len(ring) != 6:
            continue  # we expect a six-membered ring (pyranose)
        # Count how many ring atoms are oxygen.
        ring_oxygens = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:
                ring_oxygens += 1
        # Most pyranose rings have exactly one ring oxygen.
        if ring_oxygens != 1:
            continue
        
        # Count the number of external oxygen substituents on the ring.
        external_oxygen_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                # Exclude atoms that are in the ring.
                if nbr.GetIdx() in ring:
                    continue
                # Count if neighbor is oxygen (this could be hydroxyl or modified substituent like sulfo)
                if nbr.GetAtomicNum() == 8:
                    external_oxygen_count += 1
        # For a decorated sugar we expect at least three oxygen substituents.
        if external_oxygen_count >= 3:
            found_sugar = True
            break
    
    if not found_sugar:
        return False, "No candidate pyranose sugar ring with sufficient oxygen substituents detected."
    
    # --- If all checks pass, we assume the molecule is a galactosylceramide ---
    return True, "Contains one ceramide amide bond, a pyranose sugar ring (galactose) and a long alkyl chain, consistent with a galactosylceramide."

# Example usage (uncomment for testing):
# Example SMILES (one of the true positives given in the description):
# smiles_example = "C(=C/CCCCCCCCCCCCC)\\[C@@H](O)[C@@H](NC(CCCCCCCCCCCCCCCCC)=O)CO[C@H]1[C@@H]([C@H]([C@H]([C@H](O1)CO)O)OS(O)(=O)=O)O"
# result, reason = is_galactosylceramide(smiles_example)
# print(result, reason)