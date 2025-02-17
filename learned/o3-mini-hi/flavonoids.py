"""
Classifies: CHEBI:72544 flavonoids
"""
#!/usr/bin/env python3
"""
Classifies: Flavonoids (a superclass including flavonoid, isoflavonoid, chalcone, etc.)

Heuristic:
  1. Parse the SMILES string.
  2. Require that the molecule has at least 15 carbon atoms.
  3. Require that it contains at least 2 aromatic rings.
  4. Compute the Murcko scaffold (to remove appended sugars or other large substituents).
  5. Check whether the molecule (or its scaffold) contains one of several flavonoid core SMARTS patterns.
     The patterns include flavone, isoflavone, chalcone, flavanone, flavan and aurone.
  6. As an extra safeguard, if a candidate match is found then the scaffold is checked for a minimum carbon count.
  
Note: It is well known that flavonoids can be very diverse. This heuristic does not guarantee zero false positives
or negatives but may be tuned further.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_flavonoids(smiles: str):
    """
    Determines if a molecule is (putatively) a flavonoid based on its SMILES string.
    
    The classification uses the following heuristic:
      1. Validate the SMILES string.
      2. Require that the molecule has at least 15 carbon atoms.
      3. Require that it has at least 2 aromatic rings.
      4. Compute the Murcko scaffold to remove sugars/side chains.
      5. Check if either the molecule or its scaffold contains a flavonoid core as defined by several SMARTS patterns.
         (The patterns aim to capture flavone, isoflavone, chalcone, flavanone, flavan, and aurone core motifs.)
      6. For a candidate match the scaffold must also have at least 15 carbons.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a flavonoid, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the input SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count total carbon atoms in the molecule
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(carbon_atoms) < 15:
        return False, f"Only {len(carbon_atoms)} carbon atoms found; need at least 15 for a flavonoid (C15/C16 skeleton)"
    
    # Count aromatic rings.
    # We look at each ring in the ring information and require that it is at least 6-membered and completely aromatic.
    ring_info = mol.GetRingInfo()
    aromatic_ring_count = 0
    for ring in ring_info.AtomRings():
        if len(ring) >= 6 and all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_ring_count += 1
    if aromatic_ring_count < 2:
        return False, f"Only {aromatic_ring_count} aromatic ring(s) found; flavonoids usually contain at least 2 aromatic rings"
    
    # Compute the Murcko scaffold for the molecule.
    # This often removes appended sugars and other substituents so that the flavonoid aglycone is isolated.
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    except Exception:
        scaffold = None
    
    # Helper function: count carbons in a molecule.
    def count_carbons(mol_obj):
        return sum(1 for atom in mol_obj.GetAtoms() if atom.GetAtomicNum() == 6)
    
    scaffold_carbons = count_carbons(scaffold) if scaffold is not None else 0
    
    # Define an expanded list of SMARTS patterns for flavonoid cores.
    # The patterns are heuristic. Many flavonoids share a C6-C3-C6 backbone often arranged as two phenyl rings (A and B)
    # joined by a pyran (C-ring) (flavones, flavonols, flavanones) or as an open chain (chalcones).
    # Aurones (a minor subclass) are also included.
    flavonoid_patterns = [
        ("flavone",    "c1ccc2c(c1)oc(=O)c(c2)"),                   # Typical flavone / flavonol core
        ("isoflavone", "c1ccc2c(c1)oc(=O)c(c2)-c3ccccc3"),           # Isoflavone with appended benzene ring
        ("chalcone",   "c1ccc(cc1)C(=O)C=CC2=CC=CC=C2"),             # Chalcone skeleton
        ("flavanone",  "c1ccc2c(c1)C(=O)C(O)cc2"),                  # Flavanone (saturated C-ring variant)
        ("flavan",     "c1ccc2c(c1)CC(O)cc2"),                      # Generic flavan pattern (no carbonyl)
        ("aurone",     "c1ccc2c(c1)OC(=O)C=CC2")                    # Aurone pattern: benzofuranone type
    ]
    
    # Check for a flavonoid core match in either the full molecule or its scaffold.
    for name, smarts in flavonoid_patterns:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:  # skip bad patterns
            continue
        
        # Check if the pattern is found in the original molecule or if a scaffold is available, try that too.
        if mol.HasSubstructMatch(pattern) or (scaffold is not None and scaffold.HasSubstructMatch(pattern)):
            # To help rule out small non-flavonoid cores (e.g. simple coumarins), require that the scaffold itself has at least 15 carbons.
            if scaffold is None or scaffold_carbons < 15:
                # If scaffold is missing or too small, then even a match is likely spurious.
                return False, f"Matches the {name} pattern but the core (scaffold) has too few carbons ({scaffold_carbons}), likely not a full flavonoid skeleton."
            return True, f"Matches the {name} core pattern (scaffold has {scaffold_carbons} carbons)"
    
    return False, "No recognized flavonoid core pattern detected"

# Optional test examples (only when this script is run directly)
if __name__ == "__main__":
    test_examples = {
        "4',8-Dimethylgossypetin 3-glucoside": "O1C(C(O)C(O)C(O)C1OC2=C(OC3=C(C2=O)C(O)=CC(O)=C3OC)C4=CC(O)=C(OC)C=C4)CO",
        "Erycristin": "O1C2C(C3=C1C(=C(OC)C=C3)CC=C(C)C)COC4=C2C=C(C(O)=C4)CC=C(C)C",
        "Alternanthin": "O1[C@@](CC(O)[C@H](O)[C@H]1C)(C2=C(O)C3=C(OC(=CC3=O)C4=CC(OC)=C(O)C=C4)C=C2O)[H]",
        "calycosin": "COc1ccc(cc1O)-c1coc2cc(O)ccc2c1=O",
        "wogonin": "COc1c(O)cc(O)c2c1oc(cc2=O)-c1ccccc1"
    }
    
    for name, s in test_examples.items():
        result, reason = is_flavonoids(s)
        print(f"{name}: {result} ({reason})")