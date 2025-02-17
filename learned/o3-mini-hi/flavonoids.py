"""
Classifies: CHEBI:72544 flavonoids
"""
#!/usr/bin/env python3
"""
Classifies: Flavonoids (a superclass including flavonoid, isoflavonoid, chalcone, etc.)

Heuristic:
  1. Parse the SMILES string.
  2. Check that the molecule has at least 15 carbon atoms.
  3. Check that it contains at least 2 aromatic rings.
  4. Look for substructures matching one of several flavonoid core patterns.
     The list includes patterns for flavone, isoflavone, flavanone, chalcone, and a more generic flavan motif.
If these conditions are met, the molecule is classified as a putative flavonoid.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_flavonoids(smiles: str):
    """
    Determines if a molecule is a (putative) flavonoid based on its SMILES string.
    
    The classification uses the following heuristic:
      1. Validate the SMILES string.
      2. Require that the molecule has at least 15 carbon atoms.
      3. Require that the molecule has at least 2 aromatic rings.
      4. Check if the molecule contains a substructure matching one of several known flavonoid cores.
         The SMARTS patterns include flavone, chalcone, isoflavone, flavanone, and a generic flavan.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a flavonoid, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the input SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count the carbon atoms in the molecule
    c_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(c_atoms) < 15:
        return False, f"Only {len(c_atoms)} carbon atoms found; need at least 15 for a flavonoid skeleton"
    
    # Count the number of aromatic rings
    ring_info = mol.GetRingInfo()
    aromatic_ring_count = 0
    for ring in ring_info.AtomRings():
        # Check if the ring is most likely aromatic (commonly 6-membered rings)
        if len(ring) >= 6 and all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_ring_count += 1
    if aromatic_ring_count < 2:
        return False, f"Only {aromatic_ring_count} aromatic ring(s) found; flavonoids usually contain at least 2"
    
    # Define an expanded list of SMARTS patterns for flavonoid core motifs.
    # These patterns are heuristic and intended to capture many subclasses:
    flavonoid_patterns = [
        ("flavone", "c1ccc2c(c1)oc(=O)cc2"),          # Typical flavone core (C6-C3-C6 with a ketone)
        ("isoflavone", "O=c1c(O)cc2cc(-c3ccccc3)oc2c1"), # Common isoflavone pattern
        ("chalcone", "c1ccc(cc1)C(=O)C=CC2=CC=CC=C2"),  # Chalcone backbone
        ("flavanone", "c1ccc2c(c1)C(=O)C(O)cc2"),       # Flavanone variant (saturated C-ring with carbonyl)
        ("flavan", "c1ccc2c(c1)CC(O)cc2")               # A more generic flavan pattern (lacking carbonyl)
    ]
    
    # Iterate through each pattern and see if the molecule contains a corresponding substructure.
    for name, smarts in flavonoid_patterns:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue  # If the pattern fails to parse, skip it.
        if mol.HasSubstructMatch(pattern):
            return True, f"Matches the {name} core pattern"
    
    # If none of the flavonoid core patterns match, return False.
    return False, "No recognized flavonoid core pattern detected"

# Optional test examples (only run when executing this script directly)
if __name__ == "__main__":
    test_examples = {
        "4',8-Dimethylgossypetin 3-glucoside": "O1C(C(O)C(O)C(O)C1OC2=C(OC3=C(C2=O)C(O)=CC(O)=C3OC)C4=CC(O)=C(OC)C=C4)CO",
        "Erycristin": "O1C2C(C3=C1C(=C(OC)C=C3)CC=C(C)C)COC4=C2C=C(C(O)=C4)CC=C(C)C",
        "calycosin": "COc1ccc(cc1O)-c1coc2cc(O)ccc2c1=O",
        "wogonin": "COc1c(O)cc(O)c2c1oc(cc2=O)-c1ccccc1"
    }
    
    for name, s in test_examples.items():
        result, reason = is_flavonoids(s)
        print(f"{name}: {result} ({reason})")