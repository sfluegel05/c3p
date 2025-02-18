"""
Classifies: CHEBI:134251 guaiacols
"""
#!/usr/bin/env python
"""
Classifies: guaiacols
Definition: Any phenol carrying an additional methoxy substituent 
            at the ortho-position.
This module defines a function is_guaiacols that returns True with a reason 
if an aromatic ring is found where one carbon bears a hydroxyl (-OH) group and 
an ortho (adjacent in the ring) carbon carries a methoxy (-OCH3) substituent.
"""

from rdkit import Chem

def is_guaiacols(smiles: str):
    """
    Determines if a molecule is a guaiacol (a phenol with an additional methoxy 
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
        # Add explicit hydrogens so that we can properly check for OH and CH3 groups.
        mol = Chem.AddHs(mol)
        # Sanitize to ensure aromaticity is perceived.
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Error during sanitization: {e}"
        
    # Get aromatic ring information.
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        ring_set = set(ring)
        # Iterate over atoms in the ring.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # We are only interested in aromatic carbons.
            if atom.GetSymbol() != "C" or not atom.GetIsAromatic():
                continue
            
            # Check for a phenolic -OH group attached to the ring carbon.
            has_OH = False
            for nbr in atom.GetNeighbors():
                # We want substituents outside the ring.
                if nbr.GetIdx() in ring_set:
                    continue
                if nbr.GetAtomicNum() == 8:
                    # Check that this O has at least one hydrogen to be an OH.
                    for onbr in nbr.GetNeighbors():
                        if onbr.GetAtomicNum() == 1:
                            has_OH = True
                            break
                if has_OH:
                    break
            if not has_OH:
                continue  # This carbon is not a phenolic carbon.
            
            # For a confirmed phenolic carbon, look at its ortho (ring) neighbors.
            ortho_atoms = [nbr for nbr in atom.GetNeighbors() 
                           if nbr.GetIdx() in ring_set and nbr.GetSymbol() == "C" and nbr.GetIsAromatic()]
            # Now require that at least one ortho atom has a methoxy (-OCH3) substituent.
            for ortho in ortho_atoms:
                for nbr in ortho.GetNeighbors():
                    if nbr.GetIdx() in ring_set:
                        continue  # Skip atoms that are part of the same ring.
                    # We want an oxygen substituent.
                    if nbr.GetAtomicNum() != 8:
                        continue
                    # To be a methoxy group the oxygen (O) should have exactly two neighbors:
                    # one is the ring carbon (ortho) and the other should be a methyl carbon.
                    if len(nbr.GetNeighbors()) != 2:
                        continue
                    # Identify the non-ring neighbor (candidate for methyl carbon).
                    methyl_candidate = None
                    for sub in nbr.GetNeighbors():
                        if sub.GetIdx() not in ring_set and sub.GetAtomicNum() == 6:
                            methyl_candidate = sub
                    if methyl_candidate is None:
                        continue
                    # Check that the methyl candidate is sp3 hybridized.
                    if methyl_candidate.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                        continue
                    # Check that the methyl candidate has exactly three hydrogen atoms attached.
                    # (Since we added explicit hydrogens, we can count them.)
                    num_hs = sum(1 for n in methyl_candidate.GetNeighbors() if n.GetAtomicNum() == 1)
                    if num_hs == 3:
                        return True, "Found aromatic ring with a phenolic OH and an ortho methoxy group"
    return False, "No aromatic phenol with ortho methoxy substituent found"

# Example usage (uncomment to test):
# test_smiles = [
#     "COc1ccc(O)cc1",  # guaiacol itself: 2-methoxyphenol should return True
#     "c1cc(O)ccc1"     # phenol without the methoxy should return False
# ]
# for s in test_smiles:
#     result, reason = is_guaiacols(s)
#     print(f"SMILES: {s}\nClassification: {result}\nReason: {reason}\n")