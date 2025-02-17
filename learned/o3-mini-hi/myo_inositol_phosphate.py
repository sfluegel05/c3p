"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
"""
Classifies: myo-inositol phosphate – an inositol phosphate in which the inositol component 
has myo-configuration.

This improved heuristic checks for:
  1. A valid molecule with at least one six-membered ring.
  2. The ring must consist solely of carbon atoms, each with assigned chirality.
  3. Every ring carbon must have at least one substituent (non‐ring neighbor) that is an oxygen.
     Furthermore, to avoid cases with lipid or large organic substituents, the branch starting 
     at that oxygen must contain only oxygen and phosphorus atoms (plus implicit hydrogens).
  4. At least one of these oxygen branches must contain a phosphorus atom.
  
These criteria are intended to capture the myo–inositol phosphate core while filtering out larger 
substituted molecules (e.g. inositol lipids).
"""

from rdkit import Chem

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines whether a given molecule is a myo-inositol phosphate based on its SMILES string.
    
    The heuristic applied is:
      1. The molecule must be valid.
      2. There should be at least one six-membered ring in which:
           - Every ring atom is carbon with a defined chiral center.
           - Every ring carbon has at least one substituent (atom other than those in the ring)
             that is oxygen.
           - Each oxygen substituent (examined via a short breadth-first search) must not lead into 
             any atoms other than oxygen or phosphorus.
      3. At least one oxygen substituent must eventually include a phosphorus atom.
      
    Args:
        smiles (str): SMILES representation of the molecule.
    
    Returns:
        bool: True if the molecule is classified as myo-inositol phosphate, False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure stereochemistry is assigned.
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Get ring information.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No ring system found in the molecule"
    
    # Define allowed atomic numbers in a substituent branch: oxygen (8) and phosphorus (15)
    allowed_atoms = {8, 15}
    
    # Loop over rings and search for a six-membered candidate.
    for ring in rings:
        if len(ring) != 6:
            continue  # we only check six-membered rings
        
        valid_ring = True
        branch_has_phosphorus = False
        
        # Loop over atoms in the ring.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Check that the ring atom is a carbon.
            if atom.GetAtomicNum() != 6:
                valid_ring = False
                break
            # Ensure chirality is explicitly set.
            if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
                valid_ring = False
                break
            
            # Check substituents: atoms not in the ring.
            # At least one non-ring neighbor must be present.
            non_ring_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetIdx() not in ring]
            if not non_ring_neighbors:
                valid_ring = False
                break
            
            oxygen_found = False
            for nbr in non_ring_neighbors:
                # We require that the substituent is oxygen.
                if nbr.GetAtomicNum() != 8:
                    # Found a substituent that is not oxygen; disqualify this ring.
                    valid_ring = False
                    break
                
                # Do a short breadth-first search starting at this oxygen (nbr)
                # to ensure that further atoms in the substituent branch (if any)
                # are ONLY oxygen and/or phosphorus.
                visited = set()
                queue = [nbr.GetIdx()]
                branch_valid = True
                this_branch_has_P = False
                # Limit search to a small neighborhood.
                while queue:
                    curr_idx = queue.pop(0)
                    if curr_idx in visited:
                        continue
                    visited.add(curr_idx)
                    curr_atom = mol.GetAtomWithIdx(curr_idx)
                    # Check allowed atoms; note: the branch includes the starting oxygen.
                    if curr_atom.GetAtomicNum() not in allowed_atoms:
                        branch_valid = False
                        break
                    # If a phosphorus is found along the branch, note this.
                    if curr_atom.GetAtomicNum() == 15:
                        this_branch_has_P = True
                    # Traverse neighbors, but do not go back into the ring.
                    for nb in curr_atom.GetNeighbors():
                        if nb.GetIdx() in ring:
                            continue
                        if nb.GetIdx() not in visited:
                            queue.append(nb.GetIdx())
                if not branch_valid:
                    # This substituent branch contains disallowed atoms (e.g. carbon).
                    continue  # skip this substituent branch
                # If branch is valid, then we consider that this ring atom has an acceptable oxygen.
                oxygen_found = True
                if this_branch_has_P:
                    branch_has_phosphorus = True
                # We do not require every oxygen branch to include phosphorus;
                # just at least one among all ring atoms.
            if not oxygen_found:
                # If a ring carbon does not bear any oxygen substituent of acceptable type, skip the ring.
                valid_ring = False
                break
        
        # Check that the ring had at least one oxygen branch that leads to a phosphorus.
        if valid_ring and not branch_has_phosphorus:
            valid_ring = False
        
        if valid_ring:
            return True, "Found six-membered ring of chiral carbons with oxygen substituents (only O/P allowed) and at least one branch with a phosphate attachment, consistent with myo-inositol phosphate."
    
    return False, "No six-membered inositol-like ring with appropriate chiral centers, oxygen-only substituents, and phosphate attachment was found."

# Example usage (for testing purposes):
if __name__ == "__main__":
    # Test a few myo-inositol phosphate examples:
    test_smiles_examples = [
        # 1D-myo-inositol 1,2,3,4,5-pentakisphosphate
        "O[C@H]1[C@H](OP(O)(O)=O)[C@H](OP(O)(O)=O)[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@@H]1OP(O)(O)=O",
        # 1D-myo-inositol 4-phosphate
        "O[C@@H]1[C@H](O)[C@H](O)[C@@H](OP(O)(O)=O)[C@H](O)[C@H]1O",
        # A known false positive: a lipid conjugate (should be rejected)
        "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCCCCCCCCCCC"
    ]
    
    for smi in test_smiles_examples:
        result, reason = is_myo_inositol_phosphate(smi)
        print(result, reason)