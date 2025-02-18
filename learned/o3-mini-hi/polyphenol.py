"""
Classifies: CHEBI:26195 polyphenol
"""
"""
Classifies: polyphenol
Defined as 'Members of the class of phenols that contain 2 or more benzene rings 
each of which is substituted by at least one hydroxy group.'
"""
from rdkit import Chem

def is_polyphenol(smiles: str):
    """
    Determines if a molecule is a polyphenol based on its SMILES string.
    A polyphenol is defined as a molecule that contains at least 2 benzene rings
    (aromatic 6-membered rings of carbons) with at least one free hydroxy (-OH)
    substituent on each ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a polyphenol, False otherwise.
        str: Reason explaining the classification decision.
    """
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Update properties (like ring info) to be safe.
    mol.UpdatePropertyCache()
    
    # Get the ring information from the molecule.
    ring_info = mol.GetRingInfo()  # returns an object with AtomRings() method
    atom_rings = ring_info.AtomRings()  # each ring is a tuple of atom indices

    phenol_ring_count = 0  # count of benzene rings with at least one -OH substituent

    # Iterate over all rings to find benzene rings.
    for ring in atom_rings:
        # We require the ring to have exactly 6 atoms.
        if len(ring) != 6:
            continue
        
        # Check that every atom in the ring is aromatic carbon.
        is_benzene = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6 or not atom.GetIsAromatic():
                is_benzene = False
                break
        if not is_benzene:
            continue

        # For the benzene ring, check if it has at least one -OH substituent.
        has_hydroxy = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Check every neighbor of the ring atom that is not in the ring.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                # Check that the bond is a single bond (typical for -OH groups).
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                    continue
                # Identify an oxygen with at least one hydrogen attached.
                if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() >= 1:
                    has_hydroxy = True
                    break
            if has_hydroxy:
                break

        if has_hydroxy:
            phenol_ring_count += 1

    # Based on the count of benzene rings that are substituted with a hydroxy group,
    # decide on polyphenol classification.
    if phenol_ring_count < 2:
        return False, f"Only {phenol_ring_count} benzene ring(s) with a hydroxy group found"
    else:
        return True, f"Found {phenol_ring_count} benzene rings each substituted with a hydroxy group, classifying as a polyphenol"

# For testing purposes (uncomment if needed):
# if __name__ == "__main__":
#     test_smiles = "OC1C(O)c2c(O)cc(O)cc2OC1c1cc(O)c(O)c(O)c1"  # flavan-3,3',4,4',5,5',7-heptol
#     result, reason = is_polyphenol(test_smiles)
#     print("Result:", result)
#     print("Reason:", reason)