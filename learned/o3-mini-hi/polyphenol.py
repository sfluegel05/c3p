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
    A polyphenol contains at least 2 benzene rings (c1ccccc1) with at least one hydroxy group (-OH)
    attached to each ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as a polyphenol, False otherwise
        str: Reason explaining the classification decision
    """
    # Parse the SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the benzene ring pattern using SMARTS.
    # "c1ccccc1" matches a 6-membered aromatic ring of carbons.
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    matches = mol.GetSubstructMatches(benzene_pattern, useChirality=False)

    # Use a set to store unique rings (unique sorted tuples of atom indices)
    unique_rings = set()
    for match in matches:
        unique_rings.add(tuple(sorted(match)))
        
    if len(unique_rings) < 2:
        return False, "Fewer than 2 benzene rings detected"

    # For each benzene ring (as a set of atom indices), check if it has an -OH substituent.
    # An -OH substituent: an oxygen atom (atomic number 8) connected to a carbon in the ring,
    # and that oxygen should have at least one hydrogen attached (using GetTotalNumHs() > 0).
    phenol_ring_count = 0
    for ring in unique_rings:
        has_hydroxy = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Loop over neighbors that are not in the benzene ring
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() in ring:
                    continue
                # Check if the neighbor is an oxygen with at least one hydrogen attached.
                if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() >= 1:
                    has_hydroxy = True
                    break
            if has_hydroxy:
                break
        if has_hydroxy:
            phenol_ring_count += 1

    if phenol_ring_count < 2:
        return False, f"Only {phenol_ring_count} benzene ring(s) with a hydroxy group found"
    
    return True, f"Found {phenol_ring_count} benzene rings each substituted with a hydroxy group, classifying as a polyphenol"