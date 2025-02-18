"""
Classifies: CHEBI:134251 guaiacols
"""
#!/usr/bin/env python
"""
Classifies: guaiacols
Definition: Any phenol carrying an additional methoxy substituent at the ortho-position.
We require that a benzene ring (six-membered aromatic ring) has:
  - a phenolic –OH group (an oxygen attached to the ring with exactly one hydrogen) on one carbon, and
  - at least one adjacent (ortho) ring carbon carrying a methoxy (-OCH3) substituent.
  
This improved version restricts analysis to rings where every atom is an aromatic carbon.
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
        # Add explicit hydrogens so that –OH and –OCH3 details (including hydrogen counts) are explicit.
        mol = Chem.AddHs(mol)
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Error during molecule sanitization: {e}"
    
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # Iterate only over six-membered rings that are strictly benzene: all atoms are aromatic carbons.
    for ring in rings:
        if len(ring) != 6:
            continue  # we only consider six-membered rings
        # Check that every atom in the ring is a carbon and is aromatic.
        if not all(mol.GetAtomWithIdx(idx).GetSymbol() == "C" and mol.GetAtomWithIdx(idx).GetIsAromatic() 
                   for idx in ring):
            continue
        
        # For each ring atom, record if it carries a phenolic –OH or a methoxy group.
        # We will use dictionaries with keys as atom indices in the ring.
        has_OH = {}
        has_OMe = {}
        
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            oh_found = False
            ome_found = False
            # Examine each substituent (neighbor not in the ring)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                # Case 1: Look for a hydroxyl (–OH) group.
                if nbr.GetAtomicNum() == 8:  # oxygen neighbor
                    # Count hydrogens attached to this oxygen.
                    h_neighbors = [n for n in nbr.GetNeighbors() if n.GetAtomicNum() == 1]
                    # For a phenolic OH we expect exactly one hydrogen.
                    if len(h_neighbors) == 1:
                        # Also check that, aside from the ring carbon, there is no other heavy atom attached.
                        heavy_neighbors = [n for n in nbr.GetNeighbors() if n.GetAtomicNum() != 1 and n.GetIdx() not in ring]
                        if len(heavy_neighbors) == 0:
                            oh_found = True
                # Case 2: Look for a methoxy group (-OCH3).
                if nbr.GetAtomicNum() == 8:
                    # For a methoxy group we require that this oxygen has exactly 2 neighbors: 
                    # one is the ring atom and the other should be a CH3.
                    if nbr.GetDegree() == 2:
                        other = None
                        for sub in nbr.GetNeighbors():
                            if sub.GetIdx() not in ring:
                                other = sub
                        if other is not None and other.GetAtomicNum() == 6:
                            # Check that the carbon is sp3 and has exactly 3 hydrogens.
                            if other.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
                                h_on_c = [n for n in other.GetNeighbors() if n.GetAtomicNum() == 1]
                                if len(h_on_c) == 3:
                                    ome_found = True
            has_OH[idx] = oh_found
            has_OMe[idx] = ome_found
            
        # Now, for each ring atom that has an –OH, check that one of its ortho (adjacent within the ring) atoms 
        # has a methoxy substituent.
        for idx in ring:
            if not has_OH.get(idx, False):
                continue
            # Get the ring atom object.
            atom = mol.GetAtomWithIdx(idx)
            # Look for adjacent ring atoms (ortho positions) where has_OMe is True.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring and has_OMe.get(nbr.GetIdx(), False):
                    return True, "Found benzene ring with a phenolic OH and an ortho methoxy group"
                    
    return False, "No benzene ring with ortho -OH and -OCH3 substituents found"
    
# Example usage (uncomment to test):
# test_smiles = [
#     "COc1ccc(O)cc1",  # 2-methoxyphenol (guaiacol): should return True
#     "c1cc(O)ccc1"     # phenol: should return False
# ]
# for s in test_smiles:
#     result, reason = is_guaiacols(s)
#     print(f"SMILES: {s}\nClassification: {result}\nReason: {reason}\n")