"""
Classifies: CHEBI:38131 lactol
"""
"""
Classifies: lactol
Definition: Cyclic hemiacetals formed by intramolecular addition of a hydroxy group 
to a carbonyl group. In a lactol the former carbonyl carbon becomes sp³‐hybridized
and bears one –OH and one –OR substituent (with the –OR being part of a ring).
This revised version tightens the candidate search by requiring that a candidate carbon
has at least two oxygen substituents. Then, among each pair of oxygen neighbors, one 
must carry an explicit hydrogen (–OH) and the other must be “in‐ring”, and they must 
share an allowed ring (5–7 atoms). Note that in complex carbohydrate or glycoside settings
this method may still over‐/under‐classify.
"""
from rdkit import Chem

def is_lactol(smiles: str):
    """
    Determines whether the given SMILES represents a molecule featuring a lactol center.
    
    Heuristic approach:
      1. Parse the SMILES string and add explicit hydrogens.
      2. Loop over carbon atoms that are in a ring and sp3-hybridized.
      3. For each candidate carbon, require that it has 0 or 1 hydrogen(s).
      4. For that atom, gather oxygen neighbors; then, for each pair of distinct oxygen neighbors,
         require that:
           - One oxygen (OH_candidate) has at least one explicit hydrogen (–OH group);
           - The other oxygen (OR_candidate) is in a ring (serving as –OR), and
         both appear together (with the candidate carbon) in a ring of size 5–7.
      5. If a candidate passes the above tests, return True and a reason; otherwise, return False.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if a lactol center is detected, False otherwise.
        str: A textual reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that explicit –OH groups can be reliably detected.
    mol = Chem.AddHs(mol)
    
    # Get ring information (each ring as a tuple of atom indices) and define allowed ring sizes.
    ring_info = mol.GetRingInfo().AtomRings()
    valid_ring_sizes = {5, 6, 7}
    
    # Loop over candidate carbon atoms.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        if not atom.IsInRing():
            continue
        if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue
        # For a lactol carbon (hemiacetal), allow 0 or 1 attached hydrogens.
        if atom.GetTotalNumHs() not in (0, 1):
            continue
        
        # Gather oxygen neighbors.
        oxy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
        if len(oxy_neighbors) < 2:
            continue
            
        # Iterate over distinct pairs of oxygen neighbors.
        for i in range(len(oxy_neighbors)):
            for j in range(i + 1, len(oxy_neighbors)):
                oxy1 = oxy_neighbors[i]
                oxy2 = oxy_neighbors[j]
                # Determine for each oxygen whether it carries an explicit hydrogen.
                oxy1_has_H = any(nbr.GetAtomicNum() == 1 for nbr in oxy1.GetNeighbors())
                oxy2_has_H = any(nbr.GetAtomicNum() == 1 for nbr in oxy2.GetNeighbors())
                # Determine whether each oxygen is situated in a ring.
                oxy1_in_ring = oxy1.IsInRing()
                oxy2_in_ring = oxy2.IsInRing()
                
                # Assign roles: one oxygen should serve as –OH (with H) and the other as –OR (in ring).
                if oxy1_has_H and oxy2_in_ring:
                    OH_candidate = oxy1
                    OR_candidate = oxy2
                elif oxy2_has_H and oxy1_in_ring:
                    OH_candidate = oxy2
                    OR_candidate = oxy1
                else:
                    continue  # this pair fails the role criteria
                
                # Verify that the candidate carbon and the OR_candidate appear together
                # in at least one ring whose size is allowed (5–7 atoms).
                carbon_idx = atom.GetIdx()
                or_idx = OR_candidate.GetIdx()
                shares_valid_ring = False
                for ring in ring_info:
                    if carbon_idx in ring and or_idx in ring and len(ring) in valid_ring_sizes:
                        shares_valid_ring = True
                        break
                if not shares_valid_ring:
                    continue

                # Candidate passes all tests; classify this molecule as featuring a lactol center.
                return True, f"Found cyclic hemiacetal (lactol) center at carbon atom index {carbon_idx}"
    
    # No candidate lactol center found.
    return False, "No cyclic hemiacetal (lactol) center found"

# Optional simple test: if run as a script, try on a known lactol compound.
if __name__ == "__main__":
    # Example: beta-ascarylopyranose (a known lactol)
    test_smiles = "C[C@@H]1O[C@H](O)[C@H](O)C[C@H]1O"
    result, reason = is_lactol(test_smiles)
    print(result, reason)