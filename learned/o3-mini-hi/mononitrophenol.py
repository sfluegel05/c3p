"""
Classifies: CHEBI:39362 mononitrophenol
"""
"""
Classifies: mononitrophenol
Definition: A mononitrophenol is defined as a phenol (an aromatic benzene ring with a single -OH substituent)
carrying exactly one nitro group ([N+](=O)[O-]) attached to the ring.
This updated approach requires that the ring is a 6-membered (benzene) ring having exactly one –OH (or its deprotonated form)
attached and one nitro group bonded to one of its carbons.
"""

from rdkit import Chem

def is_mononitrophenol(smiles: str):
    """
    Determines if a molecule is a mononitrophenol based on its SMILES string.
    
    A mononitrophenol must:
      1. Contain exactly one nitro group ([N+](=O)[O-]).
      2. Have that nitro group attached to an aromatic (benzene) ring.
      3. The benzene ring must carry exactly one hydroxyl substituent (as -OH or deprotonated -O-)
         attached directly to a ring carbon. This is used as a proxy for a phenol.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a mononitrophenol, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Add explicit hydrogens so we can count H's on oxygen atoms.
    mol = Chem.AddHs(mol)
    
    # Define SMARTS pattern for a nitro group and search.
    nitro_smarts = "[N+](=O)[O-]"
    nitro_pattern = Chem.MolFromSmarts(nitro_smarts)
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    
    # Check that there is exactly one nitro group.
    if len(nitro_matches) != 1:
        return False, f"Found {len(nitro_matches)} nitro groups; exactly one is required"
    
    # Identify the nitrogen atom in the nitro group.
    nitro_match = nitro_matches[0]  # tuple of atom indices corresponding to the nitro group.
    nitro_N_index = None
    for idx in nitro_match:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 7:
            nitro_N_index = idx
            break
    if nitro_N_index is None:
        return False, "Error detecting the nitro nitrogen atom"
    
    # Ensure the nitro group is attached to an aromatic carbon.
    nitro_N = mol.GetAtomWithIdx(nitro_N_index)
    aromatic_carbon_idx = None
    for neighbor in nitro_N.GetNeighbors():
        if neighbor.GetAtomicNum() == 6 and neighbor.GetIsAromatic():
            aromatic_carbon_idx = neighbor.GetIdx()
            break
    if aromatic_carbon_idx is None:
        return False, "Nitro group is not attached to an aromatic carbon"
    
    # Get ring information.
    rings = mol.GetRingInfo().AtomRings()
    # Identify candidate rings that contain the aromatic carbon bonded to nitro.
    candidate_rings = [ring for ring in rings if aromatic_carbon_idx in ring]
    if not candidate_rings:
        return False, "The aromatic carbon bearing the nitro group is not part of any ring"
    
    # Check candidate rings for a benzene ring (6-membered ring with all atoms as carbon)
    # and that the ring carries exactly one hydroxyl substituent.
    for ring in candidate_rings:
        if len(ring) != 6:
            continue  # not a benzene ring
        
        # Check that all atoms in the ring are carbons.
        if not all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
            continue
        
        # Count the number of hydroxyl substituents attached to the ring.
        oh_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nb in atom.GetNeighbors():
                # Consider only neighbors not in the ring.
                if nb.GetIdx() in ring:
                    continue
                # Look only at oxygen atoms.
                if nb.GetAtomicNum() != 8:
                    continue
                # Exclude oxygens that are part of a nitro group.
                # If the neighbor oxygen is in the same nitro match, skip.
                is_nitro_oxygen = False
                for nitro_idx in nitro_match:
                    if nb.GetIdx() == nitro_idx:
                        is_nitro_oxygen = True
                        break
                if is_nitro_oxygen:
                    continue
                
                # Determine if oxygen is in a hydroxyl state.
                # Either it possesses at least one hydrogen (explicitly attached)
                # or it carries a negative formal charge (phenolate form).
                h_count = nb.GetTotalNumHs()
                if h_count > 0 or nb.GetFormalCharge() == -1:
                    oh_count += 1
        if oh_count != 1:
            # For a phenol ring there must be exactly one hydroxyl substituent.
            return False, f"Ring found (with nitro attachment) has {oh_count} hydroxyl(s) instead of exactly one"

        # If we reach here then the ring is a benzene with exactly one hydroxyl substituent.
        # Our nitro group is attached to this ring as well.
        return True, "Molecule contains exactly one nitro group attached to a benzene ring with one hydroxyl substituent"
    
    # If none of the candidate rings meets the criteria then return False.
    return False, "No benzene ring with exactly one hydroxyl substituent (phenol) attached to the nitro group was found"


# Example usage for testing:
if __name__ == "__main__":
    # A simple test with 3-nitrophenol:
    test_smiles = "Oc1cccc(c1)[N+]([O-])=O"
    result, reason = is_mononitrophenol(test_smiles)
    print("Result:", result)
    print("Reason:", reason)