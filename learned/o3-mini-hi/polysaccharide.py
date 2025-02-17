"""
Classifies: CHEBI:18154 polysaccharide
"""
"""
Classifies: Polysaccharide
Definition: A biomacromolecule consisting of large numbers of monosaccharide residues linked glycosidically.
            For our purposes, we assume a polysaccharide has (roughly) more than 10 sugar rings.
This version adjusts the previous strict criteria by:
  • Requiring the overall molecular weight to be at least 500 Da.
  • Considering only rings of size 5 or 6.
  • Counting a ring as “sugar–like” if it contains at least one oxygen, if at least half of its atoms are carbons,
    and if there is at least one –OH substituent (as identified by the SMARTS "[OX2H]") attached from a ring atom.
  • Using a relaxed threshold for hydroxyl substituents (1 as opposed to 2 on 6‐membered rings).
  • Classifying the molecule as a polysaccharide if the number of sugar–like rings is at least 11.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors

def is_polysaccharide(smiles: str):
    """
    Determines if a molecule is a polysaccharide based on its SMILES string.
    
    Improvements from the previous version:
      - Checks that the molecule overall is large (MW >= 500 Da).
      - For each ring of size 5 or 6 (typical furanose/pyranose systems), the ring is considered sugar-like if:
           • The ring contains at least 1 oxygen atom.
           • At least 50% of the atoms in the ring are carbons.
           • At least one substituent off a ring atom matches the –OH pattern ("[OX2H]").
      - The molecule is classified as a polysaccharide if the count of such sugar-like rings is at least 11.
      
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        (bool, str): Tuple of decision and reason.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # First, require that the molecule is large enough to be a polysaccharide.
    mw = Descriptors.CalcExactMolWt(mol)
    if mw < 500:
        return False, f"Molecular weight too low ({mw:.1f} Da) to be a polysaccharide."
        
    # Get ring information
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No rings detected in molecule"
    
    # Pre-compile a SMARTS pattern for a free hydroxyl group
    oh_smarts = Chem.MolFromSmarts("[OX2H]")
    
    sugar_ring_count = 0
    # Process each ring
    for ring in rings:
        ring_size = len(ring)
        if ring_size not in (5, 6):  # only look at typical sugar rings
            continue
            
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Count oxygen atoms and carbon atoms in the ring
        oxygens_in_ring = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 8)
        carbons_in_ring = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 6)
        # To be sugar-like, require at least one oxygen and at least half of the atoms are carbons.
        if oxygens_in_ring < 1 or (carbons_in_ring / ring_size) < 0.5:
            continue
        
        # Count free hydroxyl substituents (i.e. an [OX2H] attached to a ring atom, not part of the ring).
        hydroxyl_count = 0
        for atom in ring_atoms:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in ring:
                    # Check if neighbor matches an -OH using the SMARTS pattern.
                    if neighbor.HasSubstructMatch(oh_smarts):
                        # Additionally, check that the bond is single.
                        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                        if bond is not None and bond.GetBondTypeAsDouble() == 1.0:
                            hydroxyl_count += 1
        # For both 5- and 6-membered rings, relax to require at least one free hydroxyl
        if hydroxyl_count >= 1:
            sugar_ring_count += 1

    threshold = 11  # More than 10 rings is the definition.
    if sugar_ring_count >= threshold:
        return True, f"Contains {sugar_ring_count} sugar-like rings (monosaccharide residues) consistent with a polysaccharide."
    else:
        return False, f"Only {sugar_ring_count} sugar-like rings detected; need at least {threshold} for polysaccharide classification."

# Example usage:
if __name__ == "__main__":
    # Test with a typical glucose ring (a single pyranose unit)
    test_smiles = "OC1OC(O)C(O)C(O)C1O"  # one pyranose unit
    result, reason = is_polysaccharide(test_smiles)
    print("Is polysaccharide:", result)
    print("Reason:", reason)