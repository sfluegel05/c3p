"""
Classifies: CHEBI:134251 guaiacols
"""
#!/usr/bin/env python
"""
Classifies: guaiacols
Definition: Any phenol carrying an additional methoxy substituent at the ortho‐position.
A valid guaiacol has a benzene ring (six-membered aromatic ring composed solely of carbons)
that has at least one –OH group (phenolic holxyl, with one hydrogen and no extra heavy atoms, with zero formal charge)
and at least one adjacent (ortho) ring atom carrying a methoxy substituent (-OCH3; where the oxygen and its CH3 fragment have zero charge and proper connectivity).
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
    
    try:
        # Add explicit hydrogens to facilitate proper counting of H atoms.
        mol = Chem.AddHs(mol)
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Error during molecule sanitization: {e}"
    
    # Get all rings and limit ourselves to six-membered aromatic rings with only carbon atoms.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    for ring in rings:
        if len(ring) != 6:
            continue
        # Check that every atom in the ring is aromatic carbon.
        if not all(mol.GetAtomWithIdx(idx).GetSymbol() == "C" and mol.GetAtomWithIdx(idx).GetIsAromatic() 
                   for idx in ring):
            continue
        
        # Build dictionaries to mark for each ring atom:
        # has_OH: True if the atom has a substituent that is a valid phenolic –OH.
        # has_OMe: True if the atom has a substituent that is a valid methoxy group.
        has_OH = {}
        has_OMe = {}
        
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            oh_found = False
            ome_found = False
            # Look into substituents: neighbors that are not part of the ring.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                # Check for phenolic –OH:
                # The neighbor oxygen must have atomic number 8, zero formal charge,
                # exactly one hydrogen connected and no additional heavy atoms.
                if nbr.GetAtomicNum() == 8 and nbr.GetFormalCharge() == 0:
                    # Count hydrogen neighbors.
                    h_neighbors = [n for n in nbr.GetNeighbors() if n.GetAtomicNum() == 1]
                    if len(h_neighbors) == 1:
                        # Count heavy neighbors outside the ring (other than current atom).
                        heavy_neighbors = [n for n in nbr.GetNeighbors() 
                                           if n.GetAtomicNum() != 1 and n.GetIdx() not in ring and n.GetIdx() != atom.GetIdx()]
                        if len(heavy_neighbors) == 0:
                            oh_found = True
                # Check for methoxy group (-OCH3):
                # The oxygen substituent should have atomic number 8, zero formal charge,
                # exactly two neighbors (one must be the ring atom and one is the CH3 group).
                if nbr.GetAtomicNum() == 8 and nbr.GetFormalCharge() == 0:
                    if nbr.GetDegree() == 2:
                        other = None
                        for sub in nbr.GetNeighbors():
                            if sub.GetIdx() not in ring:
                                other = sub
                        if other is not None:
                            # The attached carbon should be sp3, have atomic number 6, zero formal charge,
                            # and exactly three hydrogens attached.
                            if other.GetAtomicNum() == 6 and other.GetHybridization() == Chem.rdchem.HybridizationType.SP3 and other.GetFormalCharge() == 0:
                                h_on_c = [n for n in other.GetNeighbors() if n.GetAtomicNum() == 1]
                                if len(h_on_c) == 3:
                                    ome_found = True
            has_OH[idx] = oh_found
            has_OMe[idx] = ome_found
        
        # For each ring atom with a valid –OH, check if one of its adjacent (ortho) ring atoms has a valid methoxy.
        for idx in ring:
            if not has_OH.get(idx, False):
                continue
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring and has_OMe.get(nbr.GetIdx(), False):
                    return True, "Found benzene ring with a phenolic –OH and an ortho methoxy (-OCH3) group"
                    
    return False, "No benzene ring with an ortho –OH and -OCH3 substituent pattern found"

# Example usage:
if __name__ == "__main__":
    test_smiles = [
        "COc1ccc(O)cc1",  # 2-methoxyphenol (guaiacol): should return True
        "c1cc(O)ccc1"     # phenol: should return False
    ]
    for s in test_smiles:
        result, reason = is_guaiacols(s)
        print(f"SMILES: {s}\nClassification: {result}\nReason: {reason}\n")