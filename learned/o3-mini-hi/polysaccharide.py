"""
Classifies: CHEBI:18154 polysaccharide
"""
"""
Classifies: Polysaccharide
Definition: A biomacromolecule consisting of large numbers of monosaccharide residues linked glycosidically 
            (commonly more than ten residues).
The function is_polysaccharide takes a SMILES string as input and returns a boolean along with a reason.
"""

from rdkit import Chem

def is_polysaccharide(smiles: str):
    """
    Determines if a molecule is a polysaccharide based on its SMILES string.

    For our heuristic, we assume that each monosaccharide residue appears as an aliphatic five‐membered (furanose)
    or six‐membered (pyranose) ring containing exactly one oxygen and the remainder being carbons.
    If the molecule has more than 10 such rings (i.e. at least 11), we classify it as a polysaccharide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a polysaccharide, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get information about rings in the molecule
    ring_info = mol.GetRingInfo()
    if not ring_info.IsInitialized() or not ring_info.AtomRings():
        return False, "No rings detected in molecule"

    sugar_ring_count = 0
    # iterate over all rings (each ring is a tuple of atom indices)
    for ring in ring_info.AtomRings():
        ring_size = len(ring)
        if ring_size not in (5, 6):
            continue  # only consider typical sugar ring sizes
        oxygen_count = 0
        carbon_count = 0
        # Count atoms in the ring based solely on element type
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:
                oxygen_count += 1
            elif atom.GetAtomicNum() == 6:
                # also check that carbon is sp3 (aliphatic)
                if atom.GetHybridization().name == "SP3":
                    carbon_count += 1
        # For a furanose ring, expect 1 oxygen and 4 sp3 carbons.
        # For a pyranose ring, expect 1 oxygen and 5 sp3 carbons.
        if (ring_size == 5 and oxygen_count == 1 and carbon_count == 4) or \
           (ring_size == 6 and oxygen_count == 1 and carbon_count == 5):
            sugar_ring_count += 1

    if sugar_ring_count >= 11:
        return True, f"Contains {sugar_ring_count} sugar-like rings (monosaccharide residues)"
    else:
        return False, f"Only {sugar_ring_count} sugar-like rings detected; need more than 10 for polysaccharide classification"

# Example usage:
if __name__ == "__main__":
    # Example SMILES (the below is just a placeholder small sugar; real polysaccharides are huge)
    test_smiles = "OC1OC(O)C(O)C(O)C1O"  # a typical glucose pyranose unit
    result, reason = is_polysaccharide(test_smiles)
    print("Is polysaccharide:", result)
    print("Reason:", reason)