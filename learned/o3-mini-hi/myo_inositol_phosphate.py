"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
#!/usr/bin/env python
"""
Classifies: myo-inositol phosphate – an inositol phosphate in which the inositol component 
has myo-configuration.

This heuristic checks for:
  1. A valid molecule that is neutral (formal charge = 0).
  2. At least one six‐membered ring in which every ring atom:
       • Is carbon.
       • Has an explicitly assigned chiral configuration.
       • Carries at least one non‐ring substituent.
  3. For each ring atom, at least one substituent must be an oxygen branch whose atoms 
     (found via a short breadth‐first search) are limited to oxygen and phosphorus – and 
     importantly each such atom must be uncharged.
  4. At least one oxygen branch on the ring eventually contains a phosphorus atom.
  
These criteria are intended to capture a neutral myo–inositol phosphate core (as illustrated 
by the provided examples) while filtering out overly ionized (or extended lipid‐like) species.
"""

from rdkit import Chem

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines whether a given molecule is a myo-inositol phosphate based on its SMILES string.
    
    The heuristic requires that:
      - The molecule parses correctly and carries zero net formal charge.
      - There is at least one six‐membered ring in which every atom is carbon with explicitly 
        assigned chirality.
      - Every ring carbon has at least one substituent (atom not in the ring) that is oxygen.
      - For each such oxygen substituent, a short breadth‐first search confirms that the branch 
        contains only oxygen and phosphorus atoms (each with formal charge 0) and that at least one branch 
        contains phosphorus.
    
    Args:
        smiles (str): SMILES representation of the molecule.
    
    Returns:
        (bool, str): Tuple where the boolean indicates correct classification and the string 
                     provides an explanation.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Enforce neutrality: reject molecules with non-zero overall formal charge.
    if Chem.GetFormalCharge(mol) != 0:
        return False, "Molecule carries non-zero formal charge; likely a deprotonated ionized form."
    
    # Ensure that stereochemistry is assigned.
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Get ring information from the molecule.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No ring system found in the molecule"
    
    # In substituent branches we allow only oxygen (atomic number 8) and phosphorus (15),
    # but we require that these atoms be uncharged.
    allowed_atoms = {8, 15}
    
    # Loop over all rings; look for a six-membered ring meeting our criteria.
    for ring in rings:
        if len(ring) != 6:
            continue  # only consider six-membered rings
        
        valid_ring = True
        branch_has_phosphorus = False
        
        # Process each atom in the ring:
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # The ring atom must be carbon.
            if atom.GetAtomicNum() != 6:
                valid_ring = False
                break
            # The carbon must have its chiral tag explicitly assigned.
            if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
                valid_ring = False
                break
            
            # Identify non‐ring neighbors as possible substituents.
            non_ring_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetIdx() not in ring]
            if not non_ring_neighbors:
                valid_ring = False
                break
            
            oxygen_found = False
            # Examine each substituent branch.
            for nbr in non_ring_neighbors:
                # We require that the immediate neighbor be an oxygen.
                if nbr.GetAtomicNum() != 8:
                    continue  # if not oxygen, skip it
                
                # Do a BFS on this branch (starting at the oxygen) to ensure that all atoms 
                # are either oxygen or phosphorus AND that they all carry formal charge 0.
                branch_valid = True
                queue = [nbr.GetIdx()]
                visited = set()
                branch_contains_P = False
                while queue:
                    current_idx = queue.pop(0)
                    if current_idx in visited:
                        continue
                    visited.add(current_idx)
                    current_atom = mol.GetAtomWithIdx(current_idx)
                    # Check that the atom type is allowed.
                    if current_atom.GetAtomicNum() not in allowed_atoms:
                        branch_valid = False
                        break
                    # Reject atoms that are not neutrally charged.
                    if current_atom.GetFormalCharge() != 0:
                        branch_valid = False
                        break
                    if current_atom.GetAtomicNum() == 15:
                        branch_contains_P = True
                    # Traverse neighbors that are not in the ring.
                    for nb in current_atom.GetNeighbors():
                        if nb.GetIdx() in ring:
                            continue
                        if nb.GetIdx() not in visited:
                            queue.append(nb.GetIdx())
                if not branch_valid:
                    continue  # try the next substituent branch
                # Accept this branch for the ring atom.
                oxygen_found = True
                if branch_contains_P:
                    branch_has_phosphorus = True
            if not oxygen_found:
                valid_ring = False
                break
        
        # We require at least one substituent branch on the ring to eventually contain phosphorus.
        if valid_ring and branch_has_phosphorus:
            return True, ("Found six‐membered ring of chiral carbons with substituents that are limited "
                          "to uncharged oxygen/phosphorus atoms and at least one branch bearing a phosphate, "
                          "consistent with myo-inositol phosphate.")
    
    # If no ring passes the criteria return false.
    return False, "No six‐membered inositol-like ring with appropriate chiral centers, oxygen branches, and a phosphate attachment was found."

# Example usage (for testing purposes):
if __name__ == "__main__":
    test_smiles = [
        # 1D-myo-inositol 1,2,3,4,5-pentakisphosphate (neutral)
        "O[C@H]1[C@H](OP(O)(O)=O)[C@H](OP(O)(O)=O)[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@@H]1OP(O)(O)=O",
        # 1D-myo-inositol 4-phosphate (neutral)
        "O[C@@H]1[C@H](O)[C@H](O)[C@@H](OP(O)(O)=O)[C@H](O)[C@H]1O",
        # A known false positive; note the deprotonated phosphate groups:
        "[O-]P([O-])(=O)O[C@@H]1[C@H](OP([O-])([O-])=O)[C@H](OP([O-])([O-])=O)"
        " [C@@H](OP([O-])([O-])=O)[C@H](OP([O-])([O-])=O)[C@H]1OP([O-])([O-])=O"
    ]
    for smi in test_smiles:
        result, reason = is_myo_inositol_phosphate(smi)
        print(result, reason)