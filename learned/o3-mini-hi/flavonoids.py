"""
Classifies: CHEBI:72544 flavonoids
"""
#!/usr/bin/env python3
"""
Classifies: Flavonoids (a superclass of many families such as flavonoid, isoflavonoid,
neoflavonoid, chalcones, dihydrochalcones, aurones, pterocarpan, coumestans, rotenoid,
flavonolignan, homoflavonoid and flavonoid oligomers)

Heuristic:
  - First, parse the SMILES.  
  - Check that the molecule has at least 15 carbon atoms.
  - Check that it contains at least 2 aromatic rings.
  - Look for one of several flavonoid core patterns (examples include typical flavone, isoflavone,
    or chalcone motifs).
If these conditions are met, we classify the molecule as a flavonoid.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_flavonoids(smiles: str):
    """
    Determines if a molecule is a (putative) flavonoid based on its SMILES string.
    
    The classification uses a heuristic:
      1. Validates the SMILES string.
      2. Requires at least 15 carbon atoms.
      3. Requires at least 2 aromatic rings.
      4. Checks if the molecule contains a substructure matching one of several known flavonoid
         scaffolds (e.g., flavone, chalcone or isoflavone cores).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as a flavonoid, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of carbon atoms
    c_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(c_atoms) < 15:
        return False, f"Only {len(c_atoms)} carbons found; need at least 15 for a flavonoid skeleton"

    # Check for at least 2 aromatic rings in the molecule.
    ring_info = mol.GetRingInfo()
    aromatic_ring_count = 0
    for ring in ring_info.AtomRings():
        if len(ring) >= 6:  # typical aromatic ring size
            # Check if all atoms in the ring are aromatic
            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                aromatic_ring_count += 1
    if aromatic_ring_count < 2:
        return False, f"Only {aromatic_ring_count} aromatic rings found; flavonoids usually have at least 2"

    # Define several SMARTS patterns corresponding to common flavonoid cores.
    # Note: These are heuristic and may not cover every subclass.
    flavonoid_patterns = [
        ("flavone", "O=c1cc2oc(=O)cc(c2c1)"),         # typical flavone core (C6-C3-C6)
        ("chalcone", "c1ccc(cc1)C(=O)C=CC2=CC=CC=C2"),  # chalcones have two aromatic rings linked by an α,β-unsaturated carbonyl.
        ("isoflavone", "O=c1c(O)cc2cc(-c3ccccc3)oc2c1") # a common isoflavone pattern
    ]

    for name, smarts in flavonoid_patterns:
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            continue  # skip if pattern could not be parsed
        if mol.HasSubstructMatch(patt):
            return True, f"Matches the {name} core pattern"

    # If none of our flavonoid core patterns were detected, classify as non-flavonoid.
    return False, "No recognized flavonoid core pattern detected"

# Example usage:
if __name__ == "__main__":
    # A list of example SMILES strings for testing (from the provided examples).
    test_smiles = {
        "calycosin": "COc1ccc(cc1O)-c1coc2cc(O)ccc2c1=O",
        "wogonin": "COc1c(O)cc(O)c2c1oc(cc2=O)-c1ccccc1",
        "matteucinol": "COc1ccc(cc1)[C@@H]1CC(=O)c2c(O)c(C)c(O)c(C)c2O1"
    }
    for name, s in test_smiles.items():
        result, reason = is_flavonoids(s)
        print(f"{name}: {result} ({reason})")