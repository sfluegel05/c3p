"""
Classifies: CHEBI:72544 flavonoids
"""
#!/usr/bin/env python3
"""
Classifies: Flavonoids (superclass including flavonoid, isoflavonoid, chalcone, etc.)
Improved heuristic:
  1. Validate the SMILES.
  2. Require ≥15 carbon atoms and ≥2 aromatic rings in the molecule.
  3. Compute the Murcko scaffold (to remove appended sugars/side chains) and require ≥15 carbons.
  4. Count the number of aromatic rings (using ring info) in the scaffold.
     • If there are three or more aromatic rings (often indicative of a fused ring system as in flavones, etc.)
       then try to match one of the “fused-core” SMARTS (flavone, isoflavone, flavanone, flavan, aurone).
     • If exactly two aromatic rings are present (as commonly seen in chalcones) then only the chalcone pattern is used.
  5. If a pattern is found and the scaffold’s carbon count is acceptable, the molecule is identified as a flavonoid.
Note:
  This heuristic remains approximate. Even with these changes, some false positives and negatives may occur.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_flavonoids(smiles: str):
    """
    Determines if a molecule is (putatively) a flavonoid based on its SMILES string.
    
    The classification uses the following improved heuristic:
      1. Validate the SMILES.
      2. Require that the molecule has at least 15 carbon atoms.
      3. Require that the molecule has at least 2 aromatic rings.
      4. Compute the Murcko scaffold – to focus on the core structure.
      5. Count the carbon atoms in the scaffold (require ≥15)
      6. Count the number of aromatic rings in the scaffold.
         • If the scaffold shows ≥3 aromatic rings, try matching flavone, isoflavone, flavanone, flavan, or aurone patterns.
         • If the scaffold shows exactly 2 aromatic rings, try matching a chalcone pattern.
      7. If a match is found from the respective list, return True with the pattern name;
         otherwise return False with a reason.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if classified as flavonoid, False otherwise.
        str: Description of the reasoning behind the decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count total carbons in the molecule.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 15:
        return False, f"Only {total_carbons} carbon atoms found; need ≥15 for a flavonoid skeleton."
    
    # Count aromatic rings in the molecule.
    mol_ring_info = mol.GetRingInfo()
    arom_ring_count = 0
    for ring in mol_ring_info.AtomRings():
        if len(ring) >= 6 and all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            arom_ring_count += 1
    if arom_ring_count < 2:
        return False, f"Only {arom_ring_count} aromatic ring(s) found; flavonoids usually contain ≥2 aromatic rings."
    
    # Compute the Murcko scaffold to focus on the core.
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    except Exception:
        scaffold = None
    if scaffold is None:
        return False, "Could not compute Murcko scaffold."
    
    # Count carbons in the scaffold.
    scaffold_carbons = sum(1 for atom in scaffold.GetAtoms() if atom.GetAtomicNum() == 6)
    if scaffold_carbons < 15:
        return False, f"Scaffold has only {scaffold_carbons} carbons; too few for a full flavonoid skeleton."
    
    # Count aromatic rings in the scaffold.
    scaffold_ring_info = scaffold.GetRingInfo()
    scaffold_arom_rings = []
    for ring in scaffold_ring_info.AtomRings():
        if len(ring) >= 6 and all(scaffold.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            scaffold_arom_rings.append(set(ring))
    num_scaffold_arom = len(scaffold_arom_rings)
    
    # Define SMARTS for various flavonoid cores.
    # These patterns are heuristic and are defined on the scaffold.
    # For fused systems (typically ≥3 aromatic rings) we use these patterns:
    fused_patterns = {
        "flavone":    "c1cc2oc(=O)cc(c2c1)",      # typical flavone / flavonol core
        "isoflavone": "c1ccc2c(c1)oc(=O)c(c2)",    # isoflavone-like core
        "flavanone":  "c1ccc2c(c1)C(=O)[C@@H](O)cc2",  # flavanone (with a chiral center)
        "flavan":     "c1ccc2c(c1)CC(O)cc2",       # flavan (without carbonyl)
        "aurone":     "c1ccc2c(c1)OC(=O)C=CC2"     # aurone pattern (benzofuranone)
    }
    # For open-chain chalcones (commonly with 2 nonfused aromatic rings) we use a chalcone SMARTS.
    chalcone_pattern = "c1ccc(cc1)C(=O)C=CC2=CC=CC=C2"
    
    # Based on the number of aromatic rings in the scaffold, decide which set to use.
    # (If ≥3 fused aromatic rings are present => expect an intracyclic (fused) flavonoid; if only 2,
    # then consider only chalcones.)
    target_patterns = {}
    if num_scaffold_arom >= 3:
        target_patterns = fused_patterns
    elif num_scaffold_arom == 2:
        target_patterns = {"chalcone": chalcone_pattern}
    else:
        # Fall back: no clear aromatic core in scaffold.
        return False, "Scaffold does not contain a recognizable aromatic ring system."
    
    # Now try to match one of the SMARTS patterns on the scaffold.
    for name, smarts in target_patterns.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue  # skip a bad pattern
        if scaffold.HasSubstructMatch(pattern):
            return True, f"Matches the {name} core pattern (scaffold has {scaffold_carbons} carbons, {num_scaffold_arom} aromatic rings)"
    
    return False, "No recognized flavonoid core pattern detected in the scaffold."

# Optional testing (only when this script is run directly)
if __name__ == "__main__":
    # A few test examples (names and SMILES) from the provided list.
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