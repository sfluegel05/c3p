"""
Classifies: CHEBI:26195 polyphenol
"""
"""
Classifies: polyphenol
Defined as: members of the class of phenols that contain 2 or more benzene rings,
each of which is substituted by at least one hydroxy group.
Here the approach is to first detect 6-membered rings made entirely of aromatic carbons
(i.e. benzene rings). Then for each ring we check whether at least one carbon is bonded
to an oxygen which is not involved in a carbonyl (i.e. C=O) and not linked to a sulfur (avoiding sulfate).
This heuristic is designed so that even if the oxygen is further elaborated (eg glycosylated),
it is counted while excluding O– in esters or sulfates.
"""
from rdkit import Chem

def is_polyphenol(smiles: str):
    """
    Determines if a molecule is a polyphenol (per the definition above).
    We require that the molecule contains at least 2 benzene rings (6-membered aromatic C rings)
    and that each such ring has at least one oxygen substituent that is not directly involved
    in a carbonyl or sulfate group.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is classified as a polyphenol, False otherwise.
        str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol.UpdatePropertyCache()
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Count the number of benzene rings (6-membered, all aromatic carbons) that have at least one
    # substituent oxygen which does not appear to be part of a carbonyl or sulfate.
    polyphenol_ring_count = 0
    
    # Loop over each ring detected
    for ring in atom_rings:
        # Only consider 6-membered rings
        if len(ring) != 6:
            continue
        
        # Check that every atom in the ring is a carbon and aromatic.
        is_benzene = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6 or not atom.GetIsAromatic():
                is_benzene = False
                break
        if not is_benzene:
            continue
        
        # For this benzene ring, we search for an oxygen substituent that is directly attached.
        found_valid_oxy = False
        for idx in ring:
            ring_atom = mol.GetAtomWithIdx(idx)
            for nbr in ring_atom.GetNeighbors():
                # Skip if neighbor is in the ring
                if nbr.GetIdx() in ring:
                    continue
                # We are interested in oxygen neighbors
                if nbr.GetAtomicNum() == 8:
                    # Check bonds from this oxygen to see whether it is part of an unwanted motif.
                    valid = True
                    oxy_atom = nbr
                    # rule out if the oxygen is involved in a double bond to a carbon (i.e. a carbonyl)
                    for bond in oxy_atom.GetBonds():
                        # if the bond order is DOUBLE and the other atom is carbon
                        if bond.GetBondTypeAsDouble() >= 2.0:
                            other = bond.GetOtherAtom(oxy_atom)
                            if other.GetAtomicNum() == 6:
                                valid = False
                                break
                    if not valid:
                        continue
                    # rule out oxygen that is attached to sulfur (as in sulfate groups)
                    for neighbor in oxy_atom.GetNeighbors():
                        if neighbor.GetAtomicNum() == 16:
                            valid = False
                            break
                    if valid:
                        found_valid_oxy = True
                        break
            if found_valid_oxy:
                break
        
        if found_valid_oxy:
            polyphenol_ring_count += 1

    if polyphenol_ring_count < 2:
        return False, f"Only {polyphenol_ring_count} benzene ring(s) with a valid hydroxy substituent found"
    else:
        return True, f"Found {polyphenol_ring_count} benzene rings each with a valid hydroxy substituent, classifying as a polyphenol"

# Example usage (uncomment to test locally):
# if __name__ == "__main__":
#     test_examples = [
#         # True examples:
#         "O1C2=C(OC)C(=C(C3=C(O)C=4OCOC4C(=C3C)OC)C(=C2OC1)O)C",  # Benzocamphorin E
#         "OC1C(O)c2c(O)cc(O)cc2OC1c1cc(O)c(O)c(O)c1",  # flavan-3,3',4,4',5,5',7-heptol
#         # (Add additional test SMILES here)
#         # False example:
#         "O=C(OC1=C(C(O)=C(C(=O)O)C(=C1C)C)C)C2=C(OC)C(=C(OC(=O)C3=C(O)C=C(O)C=C3C)C=C2C)C",  # Thielavin Z5
#     ]
#     for smi in test_examples:
#         result, reason = is_polyphenol(smi)
#         print("SMILES:", smi)
#         print("Result:", result)
#         print("Reason:", reason, "\n")