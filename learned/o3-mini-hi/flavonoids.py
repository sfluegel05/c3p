"""
Classifies: CHEBI:72544 flavonoids
"""
#!/usr/bin/env python3
"""
Classifies: Flavonoids (superclass including flavonoid, isoflavonoid, chalcone, etc.)
Improved heuristic:
  1. Validate the SMILES.
  2. Check that the molecule has at least 15 carbon atoms.
  3. Compute the Murcko scaffold to focus on the core. If the scaffold has fewer than 15 carbons,
     we fall back on the full molecule.
  4. Require that the core has between about 15 and 21 carbons.
  5. Count the aromatic rings (6-membered aromatic rings, determined from ring info). Flavonoid cores
     are expected to have either 2 (typically chalcones) or 3 fused aromatic rings (flavone, isoflavone, etc.).
  6. Based on that, use a set of SMARTS patterns (for fused cores or chalcone patterns) and try to match
     them on either the core (scaffold if available) or the full molecule.
  7. Return True if both the “skeleton” and a core pattern are satisfied; otherwise, return False with an explanation.
NOTE:
  This heuristic is approximate and improvements (such as more aggressive sugar removal) could be made.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_flavonoids(smiles: str):
    """
    Determines if a molecule is (putatively) a flavonoid based on its SMILES string.
    The improved heuristic:
      1. Validate the SMILES.
      2. Verify that the molecule (or its core scaffold) contains 15–21 carbon atoms.
      3. Require that the core contains either exactly 2 aromatic rings (typical of chalcones)
         or 3 aromatic rings (typical of fused flavonoid cores).
      4. Try matching a set of SMARTS patterns corresponding to known flavonoid substructures
         on the core (and if necessary, on the full molecule).
    Args:
      smiles (str): SMILES string of the molecule.
    Returns:
      bool: True if classified as a flavonoid, False otherwise.
      str: A reason detailing the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count total carbons in the full molecule.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 15:
        return False, f"Only {total_carbons} carbon atoms in the molecule; need at least 15 for a flavonoid skeleton."
    
    # Try to compute the Murcko scaffold.
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    except Exception:
        scaffold = None
    
    # Decide which to use as the "core": the scaffold if it has a sizable carbon count, or the entire molecule.
    core = None
    if scaffold:
        core_carbons = sum(1 for atom in scaffold.GetAtoms() if atom.GetAtomicNum() == 6)
        if core_carbons >= 15:
            core = scaffold
        else:
            core = mol
    else:
        core = mol
    
    # Count carbons in the core.
    core_carbons = sum(1 for atom in core.GetAtoms() if atom.GetAtomicNum() == 6)
    # We expect a flavonoid core to have between 15 and 21 carbon atoms.
    if not (15 <= core_carbons <= 21):
        return False, f"Core has {core_carbons} carbons (expected between 15 and 21 for flavonoid skeleton)."
    
    # Count aromatic rings in the core.
    core_ring_info = core.GetRingInfo()
    artmp = []
    for ring in core_ring_info.AtomRings():
        if len(ring) >= 6 and all(core.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            # Store as sorted tuple to later remove duplicates
            artmp.append(tuple(sorted(ring)))
    # Remove duplicate rings
    core_arom_rings = {r for r in artmp}
    num_core_arom = len(core_arom_rings)
    
    # For many flavonoids, the core is expected to have 2 aromatic rings (chalcones) or 3 (fused cores)
    if num_core_arom not in (2, 3):
        return False, f"Core has {num_core_arom} aromatic ring(s); expected 2 or 3 for a flavonoid core."
    
    # Define SMARTS patterns for flavonoid cores.
    # For 3-ring fused systems
    fused_patterns = {
        "flavone":    "c1cc2oc(=O)cc(c2c1)",      # typical flavone / flavonol core
        "isoflavone": "c1ccc2c(c1)oc(=O)c(c2)",    # isoflavone-like core
        "flavanone":  "c1ccc2c(c1)C(=O)[C@H](O)cc2",  # flavanone pattern (with chiral center)
        "flavan":     "c1ccc2c(c1)CC(O)cc2",       # flavan (lacks carbonyl)
        "aurone":     "c1ccc2c(c1)OC(=O)C=CC2"     # aurone pattern (benzofuranone)
    }
    # For open-chain chalcones (2 aromatic rings)
    chalcone_pattern = {"chalcone": "c1ccc(cc1)C(=O)C=CC2=CC=CC=C2"}
    
    # Based on the number of aromatic rings in the core, choose target pattern set.
    if num_core_arom == 3:
        target_patterns = fused_patterns
        target_type = "fused-core"
    elif num_core_arom == 2:
        target_patterns = chalcone_pattern
        target_type = "chalcone"
    else:
        # Should not happen because we already filtered
        return False, "Unexpected aromatic ring count in core."
    
    # We try to match the SMARTS patterns on BOTH the core and, if necessary, on the full molecule.
    # (This is to try to catch cases where the scaffold removal might have lost a key substructure.)
    matched_pattern = None
    for name, smarts in target_patterns.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue
        if core.HasSubstructMatch(pattern) or mol.HasSubstructMatch(pattern):
            matched_pattern = name
            break
    
    if matched_pattern:
        return True, f"Matches the {matched_pattern} core pattern ({target_type}) with core containing {core_carbons} carbons and {num_core_arom} aromatic rings."
    else:
        return False, "No recognized flavonoid core pattern detected in the core or full molecule."

# Optional testing block (will only run if the script is run directly)
if __name__ == "__main__":
    # A few examples drawn from the provided list.
    test_examples = {
        "4',8-Dimethylgossypetin 3-glucoside": "O1C(C(O)C(O)C(O)C1OC2=C(OC3=C(C2=O)C(O)=CC(O)=C3OC)C4=CC(O)=C(OC)C=C4)CO",
        "Erycristin": "O1C2C(C3=C1C(=C(OC)C=C3)CC=C(C)C)COC4=C2C=C(C(O)=C4)CC=C(C)C",
        "Alternanthin": "O1[C@@](CC(O)[C@H](O)[C@H]1C)(C2=C(O)C3=C(OC(=CC3=O)C4=CC(OC)=C(O)C=C4)C=C2O)[H]",
        "calycosin": "COc1ccc(cc1O)-c1coc2cc(O)ccc2c1=O",
        "wogonin": "COc1c(O)cc(O)c2c1oc(cc2=O)-c1ccccc1",
    }
    for name, s in test_examples.items():
        result, reason = is_flavonoids(s)
        print(f"{name}: {result} ({reason})")