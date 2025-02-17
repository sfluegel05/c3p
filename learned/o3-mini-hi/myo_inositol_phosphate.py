"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
"""
Classifies: myo-inositol phosphate – an inositol phosphate in which the inositol component has myo-configuration.
This heuristic checks for a six-membered saturated ring whose carbons are all individually chiral 
and each bear at least one oxygen substituent—with at least one of these oxygens being attached to a phosphorus atom.
"""

from rdkit import Chem

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines whether a given molecule is a myo-inositol phosphate based on its SMILES string.
    
    The heuristic applied is:
      1. The molecule must be valid.
      2. There should be at least one six-membered ring in which:
           - Every ring atom is carbon.
           - Every ring atom is chiral (i.e. its chiral tag is set).
           - Each ring carbon has at least one non-ring oxygen substituent.
      3. At least one of these non-ring oxygen substituents is connected to a phosphorus atom.
      
    These features are characteristic of the myo-inositol core with phosphate substituents.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a myo-inositol phosphate, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure that stereochemistry is assigned
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No ring system found in the molecule"
    
    # Loop over rings to see if any corresponds to a myo-inositol type ring
    for ring in rings:
        if len(ring) != 6:
            continue  # Only interested in six-membered rings
        
        # Check that all atoms in the ring are carbon and have defined chirality.
        all_carbons = True
        all_chiral = True
        oxygen_substituent_found = False
        phosphate_on_ring = False
        
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Check atomic number: carbon = 6
            if atom.GetAtomicNum() != 6:
                all_carbons = False
                break
            # Check chirality: must be set (i.e. not CHI_UNSPECIFIED)
            if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
                all_chiral = False
                break
            # Now: Check the neighbors that are not part of the ring.
            neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetIdx() not in ring]
            # We expect at least one oxygen (could be OH or OP... groups).
            oxygen_found_here = False
            for nbr in neighbors:
                if nbr.GetAtomicNum() == 8:
                    oxygen_found_here = True
                    # Also check if this oxygen is connected to a phosphorus (atomic number 15)
                    for subnbr in nbr.GetNeighbors():
                        if subnbr.GetAtomicNum() == 15:
                            phosphate_on_ring = True
                            break
                    # if we already found a phosphate substituent on this atom, no need to check further
                    if phosphate_on_ring:
                        break
            if oxygen_found_here:
                oxygen_substituent_found = True
            else:
                # If one ring carbon does not have an oxygen substituent, skip this ring
                oxygen_substituent_found = False
                break
        
        if not all_carbons:
            continue  # not a carbon-only ring
        if not all_chiral:
            continue  # the ring atoms do not all have defined chirality
        if not oxygen_substituent_found:
            continue  # ring atoms are missing oxygen substituents
        
        # Ensure that at least one substituent on the ring is phosphorylated.
        if not phosphate_on_ring:
            continue
        
        # If we found a valid ring, we classify it as a myo-inositol phosphate.
        return True, "Found six-membered ring with chiral carbon atoms bearing oxygen substituents and phosphate attachment, consistent with myo-inositol phosphate."
    
    return False, "No six-membered inositol-like ring with appropriate stereo centers and phosphate substituents was found."

# Example usage (for testing purposes):
if __name__ == "__main__":
    test_smiles = "O[C@H]1[C@H](O)[C@H](O)[C@@H](OP(O)(O)=O)[C@H](O)[C@@H]1O"  # 1D-myo-inositol 4-phosphate example
    result, reason = is_myo_inositol_phosphate(test_smiles)
    print(result, reason)