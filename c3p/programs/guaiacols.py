"""
Classifies: CHEBI:134251 guaiacols
"""
#!/usr/bin/env python
"""
Classifies: guaiacols
Definition: Any phenol carrying an additional methoxy substituent at the ortho‐position.
A valid guaiacol (by our definition) has a six‐membered aromatic ring (only carbon atoms) that is not fused
to a second ring and that bears at least one –OH group (free phenolic hydroxyl, with exactly one H attached)
and at least one adjacent (ortho) ring atom bearing an –OCH3 substituent (with the O attached only to the ring
and a CH3 group, where the CH3 carbon is sp3 and carries three hydrogens).
"""

from rdkit import Chem

def is_guaiacols(smiles: str):
    """
    Determines if a molecule is classified as a guaiacol (a phenol with an additional methoxy
    substituent at an ortho-position) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a guaiacol, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens and sanitize the molecule.
    try:
        mol = Chem.AddHs(mol)
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Error during molecule sanitization: {e}"
    
    # Get ring information.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # We will restrict to six-membered rings of aromatic carbon atoms.
    for ring in rings:
        if len(ring) != 6:
            continue
        # Check that every atom in the ring is an aromatic carbon.
        if not all(mol.GetAtomWithIdx(idx).GetSymbol() == "C" and mol.GetAtomWithIdx(idx).GetIsAromatic() 
                   for idx in ring):
            continue

        # To reduce false positives we require that the ring is isolated,
        # i.e. none of its atoms are shared with another six-membered aromatic ring.
        shared = False
        for other_ring in rings:
            if other_ring == ring:
                continue
            # If there is an overlap of two or more atoms, we consider the ring fused.
            if len(set(ring) & set(other_ring)) > 1:
                shared = True
                break
        if shared:
            continue

        # Dictionaries to hold if a given ring atom carries a valid free hydroxyl (-OH)
        # or a valid methoxy (-OCH3) substituent.
        has_OH = {}
        has_OMe = {}
        
        # For each atom in the ring, look at its neighbors not in the ring.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            oh_valid = False
            ome_valid = False
            
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue

                # Check for phenolic –OH:
                # The substituent oxygen must have atomic number 8, zero formal charge,
                # exactly one hydrogen attached (explicitly) and one bond (to the ring atom).
                if nbr.GetAtomicNum() == 8 and nbr.GetFormalCharge() == 0:
                    # Use GetTotalNumHs(explicitOnly=True) so that we count only added hydrogens.
                    if nbr.GetDegree() == 2 and nbr.GetTotalNumHs(explicitOnly=True) == 1:
                        # Also ensure that besides the ring atom, no other heavy atom is attached.
                        heavy_neighbors = [n for n in nbr.GetNeighbors() 
                                           if n.GetAtomicNum() != 1 and n.GetIdx() != idx]
                        if len(heavy_neighbors) == 0:
                            oh_valid = True

                # Check for methoxy group (-OCH3):
                # The oxygen must have atomic number 8, zero formal charge and exactly two bonds.
                if nbr.GetAtomicNum() == 8 and nbr.GetFormalCharge() == 0:
                    if nbr.GetDegree() == 2:
                        # One neighbor must be the ring atom; find the other.
                        other = None
                        for sub in nbr.GetNeighbors():
                            if sub.GetIdx() != idx:
                                other = sub
                        if other is not None:
                            # The attached carbon (the CH3) must be sp3, have atomic number 6,
                            # zero formal charge and exactly three hydrogens (counting explicit ones).
                            if (other.GetAtomicNum() == 6 and 
                                other.GetHybridization() == Chem.rdchem.HybridizationType.SP3 and 
                                other.GetFormalCharge() == 0):
                                if other.GetTotalNumHs(explicitOnly=True) == 3:
                                    ome_valid = True
            has_OH[idx] = oh_valid
            has_OMe[idx] = ome_valid
        
        # Now, for each ring atom that carries a free –OH,
        # check if one of its ring neighbors (ortho positions) carries a valid methoxy.
        for idx in ring:
            if not has_OH.get(idx, False):
                continue
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring and has_OMe.get(nbr.GetIdx(), False):
                    return True, ("Found isolated benzene ring with a free phenolic –OH and an "
                                  "ortho methoxy (-OCH3) substituent")
                    
    return False, "No isolated benzene ring with an ortho –OH and -OCH3 substituent pattern found"


# Example usage:
if __name__ == "__main__":
    # Two small test cases: 2-methoxyphenol (guaiacol) and phenol.
    test_smiles = [
        "COc1ccc(O)cc1",  # 2-methoxyphenol (true)
        "c1cc(O)ccc1"     # phenol (false)
    ]
    
    for s in test_smiles:
        result, reason = is_guaiacols(s)
        print(f"SMILES: {s}\nClassification: {result}\nReason: {reason}\n")