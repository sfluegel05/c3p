"""
Classifies: CHEBI:39362 mononitrophenol
"""
"""
Classifies: mononitrophenol
Definition: A nitrophenol is defined as a phenol (an aromatic ring with an –OH substituent) 
carrying a single nitro substituent ([N+](=O)[O-]) attached to the ring.
"""

from rdkit import Chem

def is_mononitrophenol(smiles: str):
    """
    Determines if a molecule is a mononitrophenol based on its SMILES string.
    A mononitrophenol should contain exactly one nitro group attached to an aromatic ring 
    that also carries a hydroxyl substituent (–OH) directly bonded to one of its aromatic atoms.
    
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

    # Define SMARTS pattern for the nitro group.
    nitro_smarts = "[N+](=O)[O-]"
    nitro_pattern = Chem.MolFromSmarts(nitro_smarts)
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    
    # Check that there is exactly one nitro group.
    if len(nitro_matches) != 1:
        return False, f"Found {len(nitro_matches)} nitro groups; exactly one is required"

    # Identify the nitrogen atom in the nitro group (atomic number 7).
    nitro_match = nitro_matches[0]  # tuple of atom indices that match the pattern.
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
    aromatic_carbon = None
    for neighbor in nitro_N.GetNeighbors():
        if neighbor.GetAtomicNum() == 6 and neighbor.GetIsAromatic():
            aromatic_carbon = neighbor.GetIdx()
            break
    if aromatic_carbon is None:
        return False, "Nitro group is not attached to an aromatic carbon"

    # Get rings in the molecule.
    rings = mol.GetRingInfo().AtomRings()
    # Identify candidate rings that include the aromatic carbon attached to the nitro group.
    candidate_rings = [ring for ring in rings if aromatic_carbon in ring]
    if not candidate_rings:
        return False, "The aromatic carbon bearing the nitro group is not part of any ring"

    # For each candidate ring, check if there is a hydroxyl (-OH) substituent.
    # A hydroxyl group is identified as an oxygen (atomic number 8) attached to a ring atom,
    # with at least one hydrogen attached (using GetTotalNumHs()).
    phenol_found = False
    for ring in candidate_rings:
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Look at neighbors that are not in the ring.
            for nb in atom.GetNeighbors():
                if nb.GetIdx() not in ring and nb.GetAtomicNum() == 8:
                    # Check that the oxygen is in a hydroxyl state by verifying it has attached hydrogens.
                    if nb.GetTotalNumHs() > 0:
                        phenol_found = True
                        break
            if phenol_found:
                break
        if phenol_found:
            break

    if not phenol_found:
        return False, "No attached hydroxyl (-OH) group detected on the aromatic ring with the nitro substituent"

    return True, "Molecule contains exactly one nitro group attached to a phenol ring"

# Example usage for testing:
if __name__ == "__main__":
    # Test with 3-nitrophenol: SMILES "Oc1cccc(c1)[N+]([O-])=O"
    test_smiles = "Oc1cccc(c1)[N+]([O-])=O"
    result, reason = is_mononitrophenol(test_smiles)
    print("Result:", result)
    print("Reason:", reason)