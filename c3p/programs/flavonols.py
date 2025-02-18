"""
Classifies: CHEBI:28802 flavonols
"""
"""
Classifies: flavonols
Definition: Any hydroxyflavone in which the ring hydrogen at position 3
of the heterocyclic C ring (of the 2-phenylchromen-4-one framework) is replaced 
by a hydroxy group (free or glycosylated). 

Our approach is two–step:
  1. Use a looser SMARTS query to detect the flavone (2-phenylchromen-4‐one) core.
     We use: "c1cc2oc(=O)c(c2c1)" – this generally matches a benzopyran-4–one.
  2. From the matched “core” we inspect the fused heterocycle (the C ring).
     In the C ring there is a ring oxygen (in the pyran) that connects two carbons.
     One of them (C2) is fused to the B ring (the extra phenyl) while the other (C3)
     should bear the –OH (or O–glycoside) group. We identify these two neighbors by 
     checking which one has an extra substituent outside the core (the B ring). Then 
     C3 is the one that lacks an extra (non–core) neighbor.
  3. Finally, we check that the candidate C3 has at least one extra substituent that 
     is an oxygen. In our heuristic, if that oxygen is bound only to a methyl group 
     (and nothing else), we consider it a methoxy and reject it.
     
Note: This is a heuristic approach and may fail in extreme cases.
"""

from rdkit import Chem

def is_flavonols(smiles: str):
    """
    Determines whether a molecule is a flavonol (3-hydroxyflavone / 3-O–glycoside).
    
    This function first confirms the presence of a flavone (2-phenylchromen-4-one) 
    core, then inspects the pyran (C) ring to see whether the carbon at position 3 
    (the one not fused to the external phenyl, i.e. the “B ring”) bears an oxygen substituent.
    
    Args:
      smiles (str): SMILES representation of the molecule.
    
    Returns:
      bool: True if the molecule is a flavonol, False otherwise.
      str: Explanation of the classification decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Use a loose SMARTS to detect the flavone core (benzopyran-4-one).
    # The pattern "c1cc2oc(=O)c(c2c1)" matches a 2-phenylchromen-4-one system.
    core_smarts = "c1cc2oc(=O)c(c2c1)"
    core_query = Chem.MolFromSmarts(core_smarts)
    if core_query is None:
        return False, "Error generating SMARTS query for flavone core"
    
    matches = mol.GetSubstructMatches(core_query)
    if not matches:
        return False, "Molecule does not contain a recognizable flavone (benzopyran-4-one) core"
    
    # For simplicity, we work on the first match.
    match = matches[0]
    core_idx_set = set(match)
    
    # From the matched core, find the heterocyclic oxygen that belongs to the pyran ring.
    # In the query, our core has two oxygen atoms: one in the ring (lowercase 'o')
    # and one as the carbonyl attached to a carbon. We choose the oxygen that is in a ring.
    ring_ox_idx = None
    for idx in match:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetSymbol().lower() == "o" and atom.IsInRing():
            # We assume the ring oxygen in the fused C ring has degree 2.
            if atom.GetDegree() == 2:
                ring_ox_idx = idx
                break
    if ring_ox_idx is None:
        return False, "Could not identify the heterocyclic (ring) oxygen of the flavone core"
    
    ring_ox_atom = mol.GetAtomWithIdx(ring_ox_idx)
    # Its neighbors (which should be carbons from the C ring) are our candidates for C2 and C3.
    nbrs = [nbr for nbr in ring_ox_atom.GetNeighbors() if nbr.GetIdx() in core_idx_set and nbr.GetAtomicNum() == 6]
    if len(nbrs) != 2:
        return False, "Unexpected number of carbon neighbors for the ring oxygen"
    
    # Heuristic to decide which neighbor is C2 vs C3:
    # The carbon fused to the external phenyl ring (B ring) is expected to have an extra neighbor 
    # that is not part of the core (i.e. part of the B ring). We loop over the two neighbors.
    candidate_C2 = None
    candidate_C3 = None
    for nbr in nbrs:
        external_neighbors = [a for a in nbr.GetNeighbors() if a.GetIdx() not in core_idx_set]
        # If a neighbor has an external connection (and the atom is aromatic) we assume it is C2.
        if external_neighbors:
            candidate_C2 = nbr
        else:
            candidate_C3 = nbr
    # If the heuristic fails (both show external connectivity or none), choose by count:
    if candidate_C2 is None and candidate_C3 is None:
        # Fall back: pick the neighbor with more external connections as C2.
        if nbrs[0].GetNumNeighbors() >= nbrs[1].GetNumNeighbors():
            candidate_C2, candidate_C3 = nbrs[0], nbrs[1]
        else:
            candidate_C2, candidate_C3 = nbrs[1], nbrs[0]
    elif candidate_C2 is None:
        candidate_C2 = nbrs[0] if nbrs[0] != candidate_C3 else nbrs[1]
    elif candidate_C3 is None:
        candidate_C3 = nbrs[0] if nbrs[0] != candidate_C2 else nbrs[1]
    
    # Now candidate_C3 should be the carbon (position 3) that must bear the hydroxy substituent.
    # Check for an oxygen substituent on candidate_C3 that is not part of the core.
    valid_3OH_found = False
    explanation_details = ""
    for subnbr in candidate_C3.GetNeighbors():
        if subnbr.GetIdx() in core_idx_set:
            continue
        if subnbr.GetAtomicNum() == 8:  # oxygen attached
            # Gather information on this oxygen's neighbors (outside the core).
            ext_neighbors = [a for a in subnbr.GetNeighbors() if a.GetIdx() not in core_idx_set]
            # If the oxygen is directly connected to a hydrogen, we assume it is a free OH.
            has_H = any(neigh.GetAtomicNum() == 1 for neigh in subnbr.GetNeighbors())
            # Otherwise if it has one connection to a carbon (and that carbon is a CH3 group),
            # we might be looking at a methoxy rather than hydroxy.
            if not has_H and len(ext_neighbors) == 1:
                carbon = ext_neighbors[0]
                # Check if the carbon is CH3 (i.e. it has three hydrogens and no other heavy atom).
                if carbon.GetAtomicNum() == 6 and carbon.GetTotalNumHs() == 3:
                    explanation_details = "The 3-O substituent appears to be methoxy (OCH3)"
                    continue  # not acceptable
            # Otherwise, assume that any oxygen substituent (or glycosylation branch) is acceptable.
            valid_3OH_found = True
            explanation_details = "Found a substituent at C3 with an oxygen (free or glycosylated)"
            break

    if not valid_3OH_found:
        return False, "Flavone core is found but the C3 carbon does not bear a valid hydroxy/glycosyl substituent"
    
    return True, f"Molecule contains a valid 3-hydroxyflavone (flavonol) core; {explanation_details}"


# Example testing code (these can be removed or modified as needed)
if __name__ == "__main__":
    # Test with one of the provided examples: tambulin
    test_smiles = "COc1ccc(cc1)-c1oc2c(OC)c(OC)cc(O)c2c(=O)c1O"
    result, reason = is_flavonols(test_smiles)
    print("SMILES:", test_smiles)
    print("Classification:", result)
    print("Reason:", reason)