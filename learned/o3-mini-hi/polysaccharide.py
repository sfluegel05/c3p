"""
Classifies: CHEBI:18154 polysaccharide
"""
"""
Classifies: Polysaccharide
Definition: A biomacromolecule consisting of large numbers of monosaccharide residues linked glycosidically 
            (commonly more than ten residues).
The function is_polysaccharide takes a SMILES string as input and returns a tuple (bool, reason).
"""

from rdkit import Chem

def is_polysaccharide(smiles: str):
    """
    Determines if a molecule is a polysaccharide based on its SMILES string.

    For our heuristic, we assume that each monosaccharide residue appears as an aliphatic five‐membered (furanose)
    or six‐membered (pyranose) ring containing exactly one oxygen and the remainder being sp3 carbons.
    If the molecule has more than 10 such rings (i.e. at least 11), we classify it as a polysaccharide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a polysaccharide, False otherwise.
        str: A reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # returns a tuple of tuples of atom indices
    if not rings:
        return False, "No rings detected in molecule"

    sugar_ring_count = 0
    # iterate over all rings
    for ring in rings:
        ring_size = len(ring)
        # only consider typical sugar ring sizes (5 or 6)
        if ring_size not in (5, 6):
            continue
        oxygen_count = 0
        carbon_count = 0
        # Count atoms by element type (we only check for oxygen and sp3 carbon)
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            atomic_num = atom.GetAtomicNum()
            if atomic_num == 8:
                oxygen_count += 1
            elif atomic_num == 6:
                if atom.GetHybridization().name == "SP3":
                    carbon_count += 1
        # For a furanose (5-membered ring): expect 1 oxygen and 4 carbons.
        # For a pyranose (6-membered ring): expect 1 oxygen and 5 carbons.
        if (ring_size == 5 and oxygen_count == 1 and carbon_count == 4) or \
           (ring_size == 6 and oxygen_count == 1 and carbon_count == 5):
            sugar_ring_count += 1

    if sugar_ring_count >= 11:
        return True, f"Contains {sugar_ring_count} sugar-like rings (monosaccharide residues) consistent with a polysaccharide."
    else:
        return False, f"Only {sugar_ring_count} sugar-like rings detected; need more than 10 for polysaccharide classification."

# Example usage:
if __name__ == "__main__":
    # This is an example; note that realistic polysaccharides are large.
    # Here we test with a single glucose pyranose unit SMILES.
    test_smiles = "OC1OC(O)C(O)C(O)C1O"  # a typical glucose unit (pyranose)
    result, reason = is_polysaccharide(test_smiles)
    print("Is polysaccharide:", result)
    print("Reason:", reason)