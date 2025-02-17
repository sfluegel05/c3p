"""
Classifies: CHEBI:26195 polyphenol
"""
"""
Classifies: polyphenol
Defined as: members of the class of phenols that contain 2 or more benzene rings,
each of which is substituted by at least one hydroxy (–OH) group.
The approach here is to add explicit hydrogens and then identify 6-membered aromatic rings (benzene rings).
For each benzene ring we check if at least one carbon in the ring is directly bonded (via a single bond)
to an oxygen that carries a hydrogen (i.e. an –OH group) and is not part of a carbonyl or sulfate.
This heuristic is designed to recognize the phenol motif.
"""
from rdkit import Chem

def is_polyphenol(smiles: str):
    """
    Determines if a molecule is a polyphenol based on the following heuristic:
    • The molecule must have at least 2 benzene rings (6-membered aromatic rings composed solely of carbons).
    • Each such ring must have at least one carbon that is substituted by an –OH group (an oxygen atom
      that is connected by a single bond and that bears at least one hydrogen).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a polyphenol, False otherwise.
        str: Explanation for the decision.
    """
    # Parse the SMILES to get a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that -OH groups are recognized properly.
    mol = Chem.AddHs(mol)
    # Update properties (e.g. aromaticity) after adding hydrogens.
    Chem.SanitizeMol(mol)
    
    # Get ring information from the molecule.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Counter for rings that appear to be benzene rings with at least one free hydroxy group.
    valid_benzene_ring_count = 0
    
    # Iterate over each ring detected by RDKit.
    for ring in atom_rings:
        # Only consider six-membered rings.
        if len(ring) != 6:
            continue
        
        # Verify that every atom in this ring is an aromatic carbon.
        is_benzene = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6 or not atom.GetIsAromatic():
                is_benzene = False
                break
        if not is_benzene:
            continue
        
        # For this aromatic ring, check if at least one ring carbon has a directly attached –OH group.
        found_hydroxy = False
        for idx in ring:
            ring_atom = mol.GetAtomWithIdx(idx)
            for nbr in ring_atom.GetNeighbors():
                # Skip atoms that are part of the ring.
                if nbr.GetIdx() in ring:
                    continue
                # Interested only in oxygen neighbors.
                if nbr.GetAtomicNum() == 8:
                    # Check that the bond between the ring atom and the oxygen is a single bond.
                    bond = mol.GetBondBetweenAtoms(ring_atom.GetIdx(), nbr.GetIdx())
                    if bond is None or bond.GetBondTypeAsDouble() != 1.0:
                        continue
                    # Check that the oxygen is actually part of an -OH group.
                    # We require that the oxygen has at least one explicit hydrogen neighbor.
                    has_H = False
                    for oxy_nbr in nbr.GetNeighbors():
                        if oxy_nbr.GetAtomicNum() == 1:
                            has_H = True
                            break
                    if has_H:
                        found_hydroxy = True
                        break
            if found_hydroxy:
                break
        
        if found_hydroxy:
            valid_benzene_ring_count += 1
    
    # A molecule is classified as a polyphenol if at least 2 benzene rings have a free -OH substituent.
    if valid_benzene_ring_count >= 2:
        return True, f"Found {valid_benzene_ring_count} benzene ring(s) each with a valid –OH substituent, classifying as a polyphenol"
    else:
        return False, f"Only {valid_benzene_ring_count} benzene ring(s) with a valid –OH substituent found"

# Example usage:
# if __name__ == "__main__":
#     test_smiles = [
#         "O([C@H]1[C@H](O)C(O[C@@H](OC2=CC=3OC=C(C(=O)C3C=C2)C4=CC=5OCOC5C=C4)C1O)CO)[C@@H]6OC([C@@H](O)[C@H](O)C6O)CO",  # Pseudobaptigenin 7-O-laminaribioside
#         "O1C2=C(OC)C(=C(C3=C(O)C=4OCOC4C(=C3C)OC)C(=C2OC1)O)C",  # Benzocamphorin E
#         "OC1C(O)c2c(O)cc(O)cc2OC1c1cc(O)c(O)c(O)c1"  # flavan-3,3',4,4',5,5',7-heptol
#     ]
#     for smi in test_smiles:
#         result, reason = is_polyphenol(smi)
#         print("SMILES:", smi)
#         print("Result:", result)
#         print("Reason:", reason)