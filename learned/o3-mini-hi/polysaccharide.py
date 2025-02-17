"""
Classifies: CHEBI:18154 polysaccharide
"""
"""
Classifies: Polysaccharide
Definition: A biomacromolecule consisting of large numbers of monosaccharide residues linked glycosidically.
            For our purposes, we assume a polysaccharide has (roughly) more than 10 sugar rings.
Improvements:
  • Uses rdMolDescriptors.CalcExactMolWt for molecular weight.
  • Only considers rings of size 5 or 6.
  • A ring is "sugar–like" if:
       – It contains at least one oxygen.
       – At least half of its atoms are carbons.
       – At least one substituent (via a single bond) attached to a ring atom is a free hydroxyl group.
         (A free hydroxyl is defined by having atomic symbol 'O' and at least one hydrogen attached.)
  • The molecule is classified as a polysaccharide if the total sugar-like rings is at least 11.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polysaccharide(smiles: str):
    """
    Determines if a molecule is a polysaccharide based on its SMILES string.
    
    Criteria:
      - The overall molecule must have a molecular weight of at least 500 Da.
      - Only rings of size 5 or 6 are considered (typical furanose or pyranose units).
      - A ring is considered "sugar–like" if:
            • It contains at least one oxygen.
            • At least half of its atoms are carbons.
            • At least one substituent (via a single bond from a ring atom) is a free hydroxyl group.
              Here, a free hydroxyl is recognized by checking that the neighbor is an oxygen atom
              (atomic number 8) and that it has one or more hydrogen atoms attached.
      - The molecule is classified as a polysaccharide if the count of sugar–like rings is at least 11.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): Tuple of decision and explanation.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Calculate the exact molecular weight using rdMolDescriptors.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 500:
        return False, f"Molecular weight too low ({mw:.1f} Da) to be a polysaccharide."
    
    # Get ring information from the molecule.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No rings detected in the molecule."
    
    sugar_ring_count = 0
    
    # Evaluate each ring.
    for ring in rings:
        ring_size = len(ring)
        # Only consider typical sugar rings: 5-membered (furanose) or 6-membered (pyranose).
        if ring_size not in (5, 6):
            continue
        
        # Get the atoms in the ring.
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Count the number of oxygen and carbon atoms in the ring.
        oxygens_in_ring = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 8)
        carbons_in_ring = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 6)
        
        # Require the ring to have at least one oxygen and at least half its atoms be carbons.
        if oxygens_in_ring < 1 or (carbons_in_ring / ring_size) < 0.5:
            continue
        
        # Count free hydroxyl substituents attached directly to the ring atoms.
        hydroxyl_count = 0
        for atom in ring_atoms:
            for neighbor in atom.GetNeighbors():
                # Only consider neighbors not in the ring (external substituents).
                if neighbor.GetIdx() in ring:
                    continue
                # Check if the neighbor atom is oxygen and if it carries at least one hydrogen.
                if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() > 0:
                    # Also ensure that the bond connecting the neighbor is a single bond.
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                    if bond is not None and bond.GetBondTypeAsDouble() == 1.0:
                        hydroxyl_count += 1
        # If at least one free hydroxyl is found on the ring, count it as sugar-like.
        if hydroxyl_count >= 1:
            sugar_ring_count += 1

    # Define the threshold—the number of sugar-like rings required.
    threshold = 11
    if sugar_ring_count >= threshold:
        return True, f"Contains {sugar_ring_count} sugar-like rings (monosaccharide residues), consistent with a polysaccharide."
    else:
        return False, f"Only {sugar_ring_count} sugar-like rings detected; need at least {threshold} for polysaccharide classification."

# Example usage:
if __name__ == "__main__":
    # Test with a single glucose unit (should be classified as NOT polysaccharide).
    test_smiles = "OC1OC(O)C(O)C(O)C1O"
    result, reason = is_polysaccharide(test_smiles)
    print("Is polysaccharide:", result)
    print("Reason:", reason)