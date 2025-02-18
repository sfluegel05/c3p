"""
Classifies: CHEBI:48024 3'-hydroxyflavanones
"""
#!/usr/bin/env python
"""
Classifies: 3'-hydroxyflavanones
Definition: Any hydroxyflavanone with a hydroxy substituent at position 3' (meta to the point of attachment on the B ring)
           of the phenyl substituent (B ring). In this version we first ensure that the molecule contains a flavanone (chromanone)
           core and then we determine the aromatic “B ring” that is attached to this core. Finally, we order the atoms in the B ring
           (which is assumed to be a benzene ring) and check that at least one of the meta positions (i.e. two atoms away on the ring)
           carries an -OH group (an oxygen with at least one hydrogen).
           
NOTE: This routine relies on a simplified SMARTS to catch the flavanone/chromanone core and so may (or may not) give the desired
      performance across all edge cases.
"""

from rdkit import Chem

def is_3__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 3'-hydroxyflavanone based on its SMILES string.
    
    Requirements:
      1. The molecule must contain a flavanone (chromanone) core. For our purpose we use a simplified SMARTS:
         "C1CC(=O)c2ccccc2O1". (This pattern is not perfect but it must be present.)
      2. At least one aromatic six-membered ring (the B ring) must be attached to any atom in the flavanone core.
      3. In that B ring, when the ring atoms are taken in order (from RDKit’s ring info) the atom that is the point
         of attachment (position 1, by our renumbering) must have at least one meta neighbor – that is, the atom at position
         (i+2 mod 6) or (i-2 mod 6) – that carries a hydroxy (-OH) group (i.e. an oxygen atom with at least one attached hydrogen).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule fits the class, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Define a simplified SMARTS for the flavanone/chromanone core.
    # This pattern represents a chroman-4-one substructure.
    core_smarts = "C1CC(=O)c2ccccc2O1"
    core_query = Chem.MolFromSmarts(core_smarts)
    if core_query is None:
        return False, "Failed to create substructure query for flavanone core"
    
    # Check if the core is present.
    core_matches = mol.GetSubstructMatches(core_query)
    if not core_matches:
        return False, "Molecule does not contain a flavanone (chromanone) core"
        
    # For each match we look for a bond from a core atom to an external aromatic atom.
    candidate_b_attachment = []  # List to hold candidate attachment atom indices on the B ring.
    core_match = core_matches[0]  # use the first occurrence of the core
    for atom_idx in core_match:
        atom = mol.GetAtomWithIdx(atom_idx)
        for nbr in atom.GetNeighbors():
            # If neighbor is not in the core and is aromatic, then consider that atom as part of the B ring.
            if nbr.GetIdx() not in core_match and nbr.GetIsAromatic():
                candidate_b_attachment.append(nbr.GetIdx())
    
    if not candidate_b_attachment:
        return False, "Could not find any aromatic substituent (B ring) attached to the flavanone core"
    
    # Get ring information.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # Now for each candidate B ring attachment we find a six-membered aromatic ring that contains it.
    for b_attach in candidate_b_attachment:
        for ring in rings:
            if len(ring) == 6 and (b_attach in ring):
                # Check all atoms in the ring are aromatic.
                if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                    continue  # not a proper benzene ring
                # Use the order provided by RDKit (ring is a tuple in ring order)
                try:
                    pos = ring.index(b_attach)
                except ValueError:
                    continue  # Should not happen
                # Determine the two meta positions: one two steps ahead and one two steps behind.
                meta_pos1 = ring[(pos + 2) % 6]
                meta_pos2 = ring[(pos - 2) % 6]
                # Check for an -OH group on these meta atoms.
                for meta_atom_idx in (meta_pos1, meta_pos2):
                    meta_atom = mol.GetAtomWithIdx(meta_atom_idx)
                    # Look through neighbors of the meta atom for an oxygen bound by a single bond
                    for nbr in meta_atom.GetNeighbors():
                        # Only consider substituents that are NOT in the B ring (to avoid counting the ring oxygen if any).
                        if nbr.GetIdx() in ring:
                            continue
                        # Check if neighbor is an oxygen.
                        if nbr.GetAtomicNum() == 8:
                            # Check if the oxygen has at least one hydrogen (typical for an -OH group).
                            # (This excludes ethers or glycosidic oxygens that lack an H.)
                            if nbr.GetTotalNumHs() >= 1:
                                # We have found an -OH group at a meta position on the B ring.
                                return True, "Molecule contains flavanone core with a -OH group at the 3' (meta) position on the B ring"
                # End of checking one candidate B ring.
    # If no candidate in any B ring carried an -OH at the meta position.
    return False, "No hydroxy (-OH) group found at the 3' (meta) position on the B ring attached to the flavanone core"

# If this script is run directly, perform a simple test.
if __name__ == "__main__":
    # Example: (2S)-3'-hydroxyflavanone (should be True)
    smiles_example = "Oc1cccc(c1)[C@@H]1CC(=O)c2ccccc2O1"
    result, reason = is_3__hydroxyflavanones(smiles_example)
    print(result, reason)