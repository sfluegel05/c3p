"""
Classifies: CHEBI:39362 mononitrophenol
"""
"""
Classifies: mononitrophenol
Definition: A mononitrophenol is defined as a phenol (a benzene ring with a single -OH substituent that is neutral)
carrying exactly one nitro group ([N+](=O)[O-]) attached to the ring.
"""

from rdkit import Chem

def is_mononitrophenol(smiles: str):
    """
    Determines if a molecule is a mononitrophenol based on its SMILES string.
    
    A mononitrophenol must:
      1. Contain exactly one nitro group ([N+](=O)[O-]).
      2. Have that nitro group attached to an aromatic benzene ring (6-membered, all carbons aromatic).
      3. This benzene ring must carry exactly one substituent hydroxyl group –OH (as a neutral oxygen with at least one hydrogen).
         (Any deprotonated oxygen, e.g. an O–, does not count.)
    
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
        
    # Add explicit hydrogens so we can properly count H's on oxygens.
    mol = Chem.AddHs(mol)
    
    # Define the SMARTS pattern for a nitro group.
    nitro_smarts = "[N+](=O)[O-]"
    nitro_pattern = Chem.MolFromSmarts(nitro_smarts)
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    
    # Check that there is exactly one nitro group present.
    if len(nitro_matches) != 1:
        return False, f"Found {len(nitro_matches)} nitro groups; exactly one is required"
    
    nitro_match = nitro_matches[0]  # tuple of atom indices in the nitro group
    
    # Identify the nitrogen atom in the nitro group.
    nitro_N = None
    for idx in nitro_match:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 7:
            nitro_N = atom
            break
    if nitro_N is None:
        return False, "Could not identify the nitro nitrogen atom"
    
    # Find an aromatic carbon neighbor attached to the nitro group.
    attached_aromatic_c = None
    for nb in nitro_N.GetNeighbors():
        if nb.GetAtomicNum() == 6 and nb.GetIsAromatic():
            attached_aromatic_c = nb
            break
    if attached_aromatic_c is None:
        return False, "Nitro group is not attached to an aromatic carbon"
    
    attached_idx = attached_aromatic_c.GetIdx()
    
    # Identify candidate rings that are benzene rings:
    # a benzene ring is defined here as a 6-membered ring where every atom is carbon and aromatic.
    rings = mol.GetRingInfo().AtomRings()
    candidate_rings = []
    for ring in rings:
        if len(ring) != 6:
            continue
        atoms_in_ring = [mol.GetAtomWithIdx(i) for i in ring]
        if not all(a.GetAtomicNum() == 6 and a.GetIsAromatic() for a in atoms_in_ring):
            continue
        # Only consider rings that include the carbon attached to the nitro group.
        if attached_idx in ring:
            candidate_rings.append(ring)
    
    if not candidate_rings:
        return False, "No benzene ring with nitro attachment found"
    
    # For each candidate benzene ring, count the hydroxyl substituents.
    # We require exactly one substituent oxygen attached to the ring that shows an O-H pattern (neutral).
    for ring in candidate_rings:
        oh_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nb in atom.GetNeighbors():
                # Consider only substituents not part of the ring.
                if nb.GetIdx() in ring:
                    continue
                # Look for oxygen atoms.
                if nb.GetAtomicNum() != 8:
                    continue
                # Skip oxygens that are part of the nitro group.
                if nb.GetIdx() in nitro_match:
                    continue
                # Check for the hydroxyl property: 
                # It must have at least one hydrogen attached and it must be neutral (formal charge 0).
                if nb.GetFormalCharge() == 0 and nb.GetTotalNumHs() >= 1:
                    oh_count += 1
                    
        if oh_count == 1:
            return True, "Molecule contains exactly one nitro group attached to a benzene ring with one neutral hydroxyl substituent"
    
    return False, "No benzene ring with exactly one neutral hydroxyl substituent (phenol) attached to the nitro group was found"


# Example usage for testing:
if __name__ == "__main__":
    # Test with 3-nitrophenol:
    test_smiles = "Oc1cccc(c1)[N+]([O-])=O"
    result, reason = is_mononitrophenol(test_smiles)
    print("Result:", result)
    print("Reason:", reason)