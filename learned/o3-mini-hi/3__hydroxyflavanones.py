"""
Classifies: CHEBI:48024 3'-hydroxyflavanones
"""
#!/usr/bin/env python
"""
Classifies: 3'-hydroxyflavanones
Definition: Any hydroxyflavanone with a hydroxy substituent at position 3' (meta to the attachment on the B ring) of the phenyl ring.
"""

from rdkit import Chem

def is_3__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 3'-hydroxyflavanone based on its SMILES string.
    
    We require:
      1. A flavanone core (a 2-phenylchroman-4-one fragment).
         Here we use the simple SMARTS "C1CC(=O)c2ccccc2O1" which captures the chromanone ring.
      2. An extra aromatic ring (“B ring”) attached to the core at the substitution carbon.
      3. That aromatic ring must have a hydroxy (-OH) substituent at the meta position
         (i.e. two bonds away from the connecting carbon of the ring).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule fits the class, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Define a SMARTS for a simplistic flavanone/chromanone core.
    # This pattern represents a 2-phenylchroman-4-one fragment.
    # The assumption is that the first atom in the match (index 0) is the chiral (or substitution) center 
    # that should be connected to the extra (B) aromatic ring.
    core_smarts = "C1CC(=O)c2ccccc2O1"
    core_query = Chem.MolFromSmarts(core_smarts)
    if core_query is None:
        return False, "Failed to create substructure query for flavanone core"
    
    core_matches = mol.GetSubstructMatches(core_query)
    if not core_matches:
        return False, "Molecule does not contain a flavanone (chromanone) core"
        
    # For this implementation, we assume that if the core is present, we can use the first match.
    core_match = core_matches[0]
    # The core SMARTS was written so that the first atom (core_match[0]) is expected to be the substitution center.
    substitution_center_idx = core_match[0]
    
    # Now, loop over neighbors of the substitution center.
    # In a flavanone, one of these neighbors (not part of the core) should be an aromatic carbon belonging to the B ring.
    atom_sub_center = mol.GetAtomWithIdx(substitution_center_idx)
    b_ring_attachment_idx = None
    for nb in atom_sub_center.GetNeighbors():
        # if this neighbor is not part of the core match then it is the attached phenyl (B) ring.
        if nb.GetIdx() not in core_match and nb.GetIsAromatic():
            b_ring_attachment_idx = nb.GetIdx()
            break
            
    if b_ring_attachment_idx is None:
        return False, "Could not find the attached aromatic (B) ring from the flavanone core"
        
    # Now we search for an aromatic ring that contains the B ring attachment atom.
    ring_info = mol.GetRingInfo()
    aromatic_rings = []
    for ring in ring_info.AtomRings():
        # We narrow to rings of size 6 (benzene rings) and all atoms aromatic.
        if len(ring) == 6:
            if b_ring_attachment_idx in ring:
                # Check that all atoms in the ring are aromatic.
                if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                    aromatic_rings.append(ring)
                    
    if not aromatic_rings:
        return False, "No aromatic six‐membered ring (B ring) found attached to the flavanone core"
    
    # Assume that the first aromatic ring that contains the attachment atom is the B ring.
    b_ring = aromatic_rings[0]
    
    # Now, we search for a hydroxy (-OH) group on the B ring 
    # at a position meta to the attachment.
    # In benzene, "meta" means the topological distance along the ring (within the ring) is 2.
    # We use RDKit’s shortest path function to measure the distance.
    found_meta_OH = False
    for atom_idx in b_ring:
        # Skip the attachment atom itself.
        if atom_idx == b_ring_attachment_idx:
            continue
        # Compute the topological distance between the attachment atom and this atom.
        path = Chem.GetShortestPath(mol, b_ring_attachment_idx, atom_idx)
        # In a benzene ring, adjacent (ortho) means path length 1, meta means 2, para means 3.
        if len(path) - 1 == 2:  # exactly meta
            # For the candidate atom, check its neighbors for an -OH substituent.
            atom = mol.GetAtomWithIdx(atom_idx)
            for nbr in atom.GetNeighbors():
                # We want an oxygen (atomic number 8) that is attached by a single bond.
                if nbr.GetAtomicNum() == 8:
                    # Using GetTotalNumHs to check if it has an attached hydrogen (typical for an -OH group).
                    # (This excludes –O–C glycosidic linkages where the oxygen would have no hydrogen.)
                    if nbr.GetTotalNumHs() >= 1:
                        found_meta_OH = True
                        break
            if found_meta_OH:
                break
                
    if not found_meta_OH:
        return False, "No hydroxy (-OH) group found at the 3' (meta) position on the B ring"
        
    # If we passed all tests:
    return True, "Molecule contains flavanone core with a hydroxy group at the 3' position on the B ring"

# For testing purposes – you could run:
if __name__ == "__main__":
    # Example: (2S)-3'-hydroxyflavanone
    smiles_example = "Oc1cccc(c1)[C@@H]1CC(=O)c2ccccc2O1"
    result, reason = is_3__hydroxyflavanones(smiles_example)
    print(result, reason)