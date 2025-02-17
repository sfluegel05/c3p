"""
Classifies: CHEBI:75874 para-terphenyl
"""
"""
Classifies: para-terphenyl derivatives
Definition: A ring assembly based on a 1,4-diphenylbenzene skeleton and its substituted derivatives thereof.
A true para-terphenyl contains a central benzene ring (non-fused, all six atoms are carbon and aromatic)
with benzene rings attached via a single (non-fused) bond at the para positions (i.e. separated by three atoms).
"""

from rdkit import Chem

def is_para_terphenyl(smiles: str):
    """
    Determines if a molecule is a para-terphenyl derivative based on its SMILES string.
    
    The algorithm does the following:
     1) Uses a SMARTS pattern ("c1ccccc1") to find all benzene rings.
     2) Filters for central ring candidates that are pure carbon.
     3) For each candidate ring, it determines an ordering (cyclic order) of its atoms.
     4) For each atom of the candidate, it looks for bonds to an external atom that belongs to
        a benzene ring (again using the benzene SMARTS) such that the external benzene ring shares
        exactly one atom with the candidate ring (to ensure non-fused attachment).
     5) It then checks whether two such bonds occur on the candidate ring in para positions 
        (i.e. opposite, 3 positions apart in a 6-membered cycle).
        
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule displays a para-terphenyl motif, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Error sanitizing molecule: {e}"
    
    # Define benzene SMARTS: a six-membered aromatic ring.
    benzene_smarts = Chem.MolFromSmarts("c1ccccc1")
    benzene_matches = mol.GetSubstructMatches(benzene_smarts)
    if not benzene_matches:
        return False, "No benzene rings found in the molecule"
    
    # Make a unique set of benzene rings (each as a frozenset of atom indices)
    benzene_rings = []
    for match in benzene_matches:
        fs = frozenset(match)
        if fs not in benzene_rings:
            benzene_rings.append(fs)
    
    # Helper: order the atoms of a six-membered ring (assumes it is a proper cycle)
    def get_cycle_ordering(mol, ring_set):
        # Convert ring_set to list so we can choose a starting atom.
        ring_atoms = set(ring_set)
        ordering = []
        # Choose an arbitrary starting atom.
        start = next(iter(ring_atoms))
        ordering.append(start)
        prev = None
        current = start
        while len(ordering) < len(ring_set):
            current_atom = mol.GetAtomWithIdx(current)
            found = False
            for nbr in current_atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx in ring_atoms and nbr_idx != prev:
                    ordering.append(nbr_idx)
                    prev = current
                    current = nbr_idx
                    found = True
                    break
            if not found:
                break  # Incomplete ordering
        return ordering
    
    # The central ring of a para-terphenyl should be a benzene ring with no heteroatoms.
    central_candidates = []
    for ring in benzene_rings:
        # Check that all atoms in the ring are carbon:
        if all(mol.GetAtomWithIdx(idx).GetSymbol() == "C" for idx in ring):
            if len(ring) == 6:  # Only proper six-membered rings are considered.
                central_candidates.append(ring)
    
    if not central_candidates:
        return False, "No pure carbon six-membered (benzene) rings found suitable as central rings"
    
    # Now check each central candidate for two para-attached benzene rings.
    for candidate in central_candidates:
        ordering = get_cycle_ordering(mol, candidate)
        if len(ordering) != 6:
            continue  # not a proper cycle; skip
        substituent_positions = []
        # For each atom (by position) in the candidate ring, look for a substituent benzene ring.
        for pos, atom_idx in enumerate(ordering):
            atom = mol.GetAtomWithIdx(atom_idx)
            found_substituent = False
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                # Only consider atoms not in candidate ring (i.e. external connection)
                if nbr_idx in candidate:
                    continue
                # Now look for any benzene ring (from our list) that includes this neighbor.
                # To be a non-fused substituent, the overlap between the candidate ring and the external benzene must be exactly one atom.
                for ext_ring in benzene_rings:
                    if nbr_idx in ext_ring and len(candidate.intersection(ext_ring)) == 1:
                        substituent_positions.append(pos)
                        found_substituent = True
                        break
                if found_substituent:
                    break
        # Must have at least two substituents on the central ring.
        if len(substituent_positions) < 2:
            continue
        # Check if any two substituents are para (difference of 3 positions in a 6-membered ring)
        n = 6
        for i in range(len(substituent_positions)):
            for j in range(i+1, len(substituent_positions)):
                diff = abs(substituent_positions[i] - substituent_positions[j])
                circular_diff = min(diff, n - diff)
                if circular_diff == 3:
                    return True, ("Found para-terphenyl skeleton: central benzene ring with two benzene substituents in para positions "
                                  "attached via single (non-fused) bonds.")
    return False, "Could not find a central benzene ring with two para-substituted benzene rings (non-fused) required for para-terphenyl"

# Example usage:
# Uncomment the following lines to test with provided SMILES examples:
# test_smiles = [
#    "O=C(OC1=C(OC(=O)C)C(=C(O)C(=C1C2=CC(O)=C(O)C=C2)O)C3=CC=C(O)C=C3)C",  # 2',3'-diacetoxy-3,4,5',6',4''-pentahydroxy-p-terphenyl
#    "O=C(OC1=C(O)C(=C(OC(=O)C[C@@H](OC(=O)C)C)C(=C1C2=CC=C(O)C=C2)O)C3=CC=C(O)C=C3)C[C@@H](OC(=O)C)C",  # Curtisian M
#    "COc1cc(-c2ccccc2)c(OC)c(O)c1-c1ccc(O)c(O)c1",  # A simplified p-terphenyl derivative
# ]
# for s in test_smiles:
#    result, reason = is_para_terphenyl(s)
#    print(f"SMILES: {s}\nResult: {result}, Reason: {reason}\n")