"""
Classifies: CHEBI:18154 polysaccharide
"""
"""
Classifies: Polysaccharide
Definition: A biomacromolecule consisting of large numbers of monosaccharide residues linked glycosidically.
            For our purposes, we assume a polysaccharide has (roughly) more than 10 sugar rings.
Improvements:
  • Uses rdMolDescriptors.CalcExactMolWt instead of Descriptors.CalcExactMolWt.
  • Checks that the overall molecular weight is at least 500 Da.
  • Reviews each ring (of size 5 or 6) to determine if it is “sugar–like” based on:
       – The ring contains at least one oxygen.
       – At least half its atoms are carbons.
       – At least one substituent attached to a ring atom matches a free hydroxyl group (“[OX2H]”).
  • Classifies as a polysaccharide if the count of such sugar–like rings is at least 11.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polysaccharide(smiles: str):
    """
    Determines if a molecule is a polysaccharide based on its SMILES string.
    
    Criteria:
      - The overall molecule must have a molecular weight of at least 500 Da.
      - Only rings of size 5 or 6 are considered (typical furanose or pyranose units).
      - A ring is "sugar–like" if:
            • It contains at least one oxygen.
            • At least half of its atoms are carbons.
            • At least one external substituent attached (via a single bond) to a ring atom matches the free –OH pattern.
      - The molecule is classified as a polysaccharide if the number of sugar-like rings is at least 11.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): Tuple of decision and explanation.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Use rdMolDescriptors to calculate the exact molecular weight.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 500:
        return False, f"Molecular weight too low ({mw:.1f} Da) to be a polysaccharide."
    
    # Fetch ring information. If no rings exist, exit.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No rings detected in molecule"
    
    # Pre-compile pattern for a free hydroxyl group.
    oh_smarts = Chem.MolFromSmarts("[OX2H]")
    sugar_ring_count = 0
    
    # Evaluate each ring.
    for ring in rings:
        ring_size = len(ring)
        # Only consider 5-membered and 6-membered rings.
        if ring_size not in (5, 6):
            continue
        
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Count oxygen and carbon atoms in the ring.
        oxygens_in_ring = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 8)
        carbons_in_ring = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 6)
        
        # Require the ring to have at least one oxygen and at least half its atoms are carbon.
        if oxygens_in_ring < 1 or (carbons_in_ring / ring_size) < 0.5:
            continue
        
        # Count free hydroxyl substituents on the ring atoms.
        hydroxyl_count = 0
        for atom in ring_atoms:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in ring:
                    # Check if the neighbor structure corresponds to a free hydroxyl.
                    if neighbor.HasSubstructMatch(oh_smarts):
                        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                        if bond is not None and bond.GetBondTypeAsDouble() == 1.0:
                            hydroxyl_count += 1
        # If at least one free hydroxyl is found, count this ring as sugar-like.
        if hydroxyl_count >= 1:
            sugar_ring_count += 1

    # For our purposes, a molecule must have at least 11 sugar-like rings.
    threshold = 11
    if sugar_ring_count >= threshold:
        return True, f"Contains {sugar_ring_count} sugar-like rings (monosaccharide residues) consistent with a polysaccharide."
    else:
        return False, f"Only {sugar_ring_count} sugar-like rings detected; need at least {threshold} for polysaccharide classification."

# Example usage:
if __name__ == "__main__":
    # Test with a simple glucose ring (not a polysaccharide).
    test_smiles = "OC1OC(O)C(O)C(O)C1O"  # A single pyranose unit.
    result, reason = is_polysaccharide(test_smiles)
    print("Is polysaccharide:", result)
    print("Reason:", reason)