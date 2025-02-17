"""
Classifies: CHEBI:134251 guaiacols
"""
#!/usr/bin/env python
"""
Classifies: guaiacols
Definition: Any phenol carrying an additional methoxy substituent at the ortho‐position.
"""
from rdkit import Chem

def is_guaiacols(smiles: str):
    """
    Determines if a molecule is a guaiacol (any phenol having an additional methoxy group at an ortho-position)
    based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a guaiacol, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Ensure that aromaticity is perceived (should be by default)
    Chem.SanitizeMol(mol)
    
    ring_info = mol.GetRingInfo()
    # Iterate over all rings in the molecule.
    for ring in ring_info.AtomRings():
        # Convert ring indices to a set for faster lookup.
        ring_set = set(ring)
        # Look in the ring for an aromatic carbon that bears a hydroxyl group.
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            # We expect phenol to be on an aromatic carbon.
            if atom.GetSymbol() != "C" or not atom.GetIsAromatic():
                continue

            # Check all neighbors of the aromatic carbon that are not in the ring.
            hydroxyl_found = False
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring_set:
                    continue  # will check substituents on the ring later for ortho methoxy.
                # Check for hydroxyl group: oxygen that carries at least one hydrogen.
                if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() > 0:
                    hydroxyl_found = True
                    break
            if not hydroxyl_found:
                continue  # Not a phenol carbon, move to next atom.
                
            # At this point, we have an aromatic carbon bearing an -OH (phenolic group).
            # Now find the ortho positions (neighbors within the same ring).
            ortho_atoms = []
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring_set and nbr.GetSymbol() == "C":
                    ortho_atoms.append(nbr)
            
            # For each ortho atom, check its non-ring substituents for a methoxy group.
            for ortho in ortho_atoms:
                for sub in ortho.GetNeighbors():
                    if sub.GetIdx() in ring_set:
                        continue  # Skip atoms in the ring.
                    # Check if the substituent is oxygen.
                    if sub.GetAtomicNum() != 8:
                        continue
                    # To be a methoxy group, the oxygen should be bonded to exactly two atoms:
                    # the ring carbon and a methyl carbon.
                    if len(sub.GetNeighbors()) != 2:
                        continue
                    # Identify the non-ring neighbor attached to this oxygen.
                    methoxy_candidate = None
                    for x in sub.GetNeighbors():
                        if x.GetIdx() not in ring_set and x.GetAtomicNum() == 6:
                            methoxy_candidate = x
                    if methoxy_candidate is None:
                        continue
                    # Check if the candidate carbon is a methyl group (has at least three hydrogens).
                    if methoxy_candidate.GetTotalNumHs() >= 3:
                        return True, "Found aromatic ring with a phenolic OH and an ortho methoxy group"
                    
    # If we exit the loop, no ring satisfied the guaiacol definition.
    return False, "No aromatic phenol with ortho methoxy substituent found"

# Example usage (uncomment the following lines to test):
# smiles_list = [
#     "COc1ccc(O)cc1",  # guaiacol itself: 2-methoxyphenol should return True
#     "c1cc(O)ccc1"     # phenol without the methoxy should return False
# ]
# for s in smiles_list:
#     result, reason = is_guaiacols(s)
#     print(f"SMILES: {s}\nClassification: {result}\nReason: {reason}\n")