"""
Classifies: CHEBI:134251 guaiacols
"""
#!/usr/bin/env python
"""
Classifies: guaiacols
Definition: Any phenol carrying an additional methoxy substituent at the ortho-position.
We require that a benzene ring (six-membered aromatic ring) has:
  - a phenolic –OH group (oxygen attached to the ring with exactly one hydrogen) on one carbon, and
  - at least one adjacent (ortho) carbon on the same ring that carries a methoxy (-OCH3) substituent,
    as defined by an oxygen connected to the ring and a methyl group (an sp³ carbon with exactly three hydrogens).
"""

from rdkit import Chem

def is_guaiacols(smiles: str):
    """
    Determines if a molecule is a guaiacol (a phenol with an additional methoxy substituent 
    at an ortho-position) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a guaiacol, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    try:
        # Add explicit hydrogens so that O–H and –OCH3 details are clear.
        mol = Chem.AddHs(mol)
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Error during sanitization: {e}"
    
    # Get ring information.
    ring_info = mol.GetRingInfo()
    
    # Iterate over all rings.
    for ring in ring_info.AtomRings():
        # Restrict search to benzene-like rings (6 members, all aromatic carbons).
        if len(ring) != 6:
            continue  # skip non-benzene rings
        # Verify that each atom in this ring is a carbon.
        if any(mol.GetAtomWithIdx(idx).GetSymbol() != "C" for idx in ring):
            continue
        
        # Create dictionaries to store which ring atoms have an OH or a methoxy group.
        has_OH = {}    # key: atom index in ring -> bool
        has_OMe = {}   # key: atom index in ring -> bool
        
        # Inspect each ring atom for substituents.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # We expect aromatic carbon already by construction.
            # Check all neighbors that are NOT part of the ring.
            oh_found = False
            ome_found = False
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue  # skip atoms in the ring
                # Case 1. Check if neighbor is –OH:
                if nbr.GetAtomicNum() == 8:
                    # Get the hydrogens attached to O.
                    h_count = sum(1 for n in nbr.GetNeighbors() if n.GetAtomicNum() == 1)
                    # For a phenolic OH we expect exactly one hydrogen.
                    if h_count == 1:
                        oh_found = True
                # Case 2. Check if neighbor is part of a methoxy group (-OCH3):
                # Criteria: oxygen atom connected directly to the ring with exactly two neighbors,
                # one in the ring (already the case) and one which is a methyl carbon.
                if nbr.GetAtomicNum() == 8:
                    if nbr.GetDegree() == 2:
                        # Identify the neighbor in addition to the ring atom.
                        other = None
                        for sub in nbr.GetNeighbors():
                            if sub.GetIdx() not in ring:
                                other = sub
                        if other is not None and other.GetAtomicNum() == 6:
                            # Check that the carbon is sp3 (methyl) and has exactly 3 hydrogens.
                            if other.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
                                h_on_c = sum(1 for n in other.GetNeighbors() if n.GetAtomicNum() == 1)
                                # Additionally, ensure that the carbon is connected only to the oxygen (and hydrogens).
                                if h_on_c == 3:
                                    ome_found = True
            has_OH[idx] = oh_found
            has_OMe[idx] = ome_found

        # Now, loop over ring atoms: if an atom has an OH, check its ortho (adjacent in ring) atoms for a methoxy.
        # For benzene, the ortho atoms are those directly bonded to this atom.
        # We check using the molecule's bond structure.
        for idx in ring:
            if not has_OH.get(idx, False):
                continue
            # Get neighbors within the ring (these are the ortho positions).
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    if has_OMe.get(nbr.GetIdx(), False):
                        return True, "Found benzene ring with a phenolic OH and an ortho methoxy group"
    return False, "No benzene ring with ortho -OH and -OCH3 substituents found"

# Example usage (uncomment to test):
# test_smiles = [
#     "COc1ccc(O)cc1",  # guaiacol itself: 2-methoxyphenol should return True
#     "c1cc(O)ccc1"     # phenol without the methoxy should return False
# ]
# for s in test_smiles:
#     result, reason = is_guaiacols(s)
#     print(f"SMILES: {s}\nClassification: {result}\nReason: {reason}\n")