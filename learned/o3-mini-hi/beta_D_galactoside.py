"""
Classifies: CHEBI:28034 beta-D-galactoside
"""
"""
Classifies: beta-D-galactoside

Definition: Any D-galactoside having beta-configuration at its anomeric centre.
A beta-D-galactoside here is defined as a molecule that contains a complete beta-D-galactopyranoside fragment.
This implementation looks for a six-membered ring (pyranose) that has exactly one ring oxygen and five ring carbons.
Within the ring, we search for a carbon with beta configuration (C@@H) that is bonded to an oxygen substituent 
(excluding ring connections) and we also require the presence of a characteristic exocyclic –CH2OH group attached 
to one of the ring carbons (expected for D-sugars).
"""

from rdkit import Chem

def is_beta_D_galactoside(smiles: str):
    """
    Determines if a molecule is a beta-D-galactoside based on its SMILES string.

    The implemented strategy consists of:
      1. Parsing the SMILES and ensuring stereochemistry is assigned.
      2. Enumerating all rings of size 6.
      3. For each 6-membered ring, ensuring the ring contains exactly one oxygen (a hallmark of a pyranose)
         and five carbons.
      4. Searching in that ring for a carbon that is chiral and has a beta configuration 
         (i.e. its chiral tag equals CHI_TETRAHEDRAL_CCW) and that is bound to an external oxygen (the glycosidic oxygen).
      5. Verifying that one (typically at the C5-equivalent position) carries a CH2OH group; this is done by requiring
         that at least one ring carbon (outside the candidate anomeric carbon) has a neighbor that is a –CH2OH fragment.
         Here we approximate a CH2OH group by checking for a carbon (atomic number 6) attached to at least one oxygen
         (atomic number 8) and not being part of the ring.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule appears to contain a beta-D-galactoside sugar moiety, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Assign stereochemistry so that chiral tags are set.
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # Loop over all rings looking for a candidate pyranose ring.
    for ring in rings:
        if len(ring) != 6:
            continue  # Only interested in 6-membered rings
        
        # Count the number of oxygen atoms in the ring.
        oxy_in_ring = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8]
        if len(oxy_in_ring) != 1:
            continue  # A pyranose ring should contain exactly one ring oxygen
        
        # Count the number of carbon atoms in the ring.
        carbons_in_ring = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6]
        if len(carbons_in_ring) != 5:
            continue  # Expect 5 carbons in a pyranose ring
        
        # Look for a candidate anomeric carbon in beta configuration.
        # Candidate anomeric carbons are ring carbons (atomic num 6) that:
        #   a) are chiral and have configuration beta (RDKit uses CHI_TETRAHEDRAL_CCW for [C@@H]),
        #   b) are attached to at least one oxygen that is not part of the ring.
        candidate_anomeric = None
        for idx in carbons_in_ring:
            atom = mol.GetAtomWithIdx(idx)
            if not atom.HasProp('_ChiralityPossible'):
                # If chirality was never set, skip.
                continue
            if atom.GetChiralTag() != Chem.CHI_TETRAHEDRAL_CCW:
                continue  # We require the beta configuration ([C@@H])
            # Check neighbors for exocyclic oxygen (not in the ring)
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8 and (nbr.GetIdx() not in ring):
                    candidate_anomeric = idx
                    break
            if candidate_anomeric is not None:
                break
        
        if candidate_anomeric is None:
            continue  # No anomeric carbon candidate found in this ring
        
        # Now check if one of the other ring carbons carries an exocyclic CH2OH group.
        # We look for a non-ring carbon neighbor (of a ring carbon) that is a carbon (atomic number 6) and is attached to at least one oxygen.
        has_CH2OH = False
        for idx in carbons_in_ring:
            # We skip the candidate anomeric carbon (which is the glycosidic linkage site)
            if idx == candidate_anomeric:
                continue
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() == 6:
                    # Check if this neighbor (potential CH2 group) has an oxygen neighbor
                    # (this is a loose check for CH2OH)
                    for subnbr in nbr.GetNeighbors():
                        if subnbr.GetAtomicNum() == 8:
                            has_CH2OH = True
                            break
                if has_CH2OH:
                    break
            if has_CH2OH:
                break
        
        if not has_CH2OH:
            continue  # The ring does not appear to have the required exocyclic CH2OH group
        
        # If we found a ring that fulfills all the requirements, we call the molecule a beta-D-galactoside.
        return True, "Contains a six-membered pyranose ring with beta-configured anomeric carbon and exocyclic CH2OH group (consistent with a beta-D-galactoside)"
    
    # If none of the rings qualifies, return a failure.
    return False, "No complete beta-D-galactoside sugar fragment found"

# Example usage (for testing):
if __name__ == "__main__":
    # One of the provided examples: methyl beta-D-galactoside should return True.
    test_smiles = "CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O"
    result, reason = is_beta_D_galactoside(test_smiles)
    print("SMILES:", test_smiles)
    print("Result:", result)
    print("Reason:", reason)