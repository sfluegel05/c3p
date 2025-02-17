"""
Classifies: CHEBI:75874 para-terphenyl
"""
"""
Classifies: para-terphenyl derivatives
Definition: A ring assembly based on a 1,4-diphenylbenzene skeleton and its substituted derivatives thereof.
A true para-terphenyl contains a central benzene ring (non-fused) that is para-substituted 
(i.e. substituents at positions 1 and 4; in a hexagon these positions are 3 bonds apart along the ring)
with benzene rings attached via a single (non-fused) bond.
"""

from rdkit import Chem

def is_para_terphenyl(smiles: str):
    """
    Determines if a molecule is a para-terphenyl derivative based on its SMILES string.
    
    The algorithm first finds all benzene rings (six-membered aromatic carbocyles) in the molecule.
    For each benzene ring candidate, it reconstructs the cyclic order of its atoms.
    Then it checks whether at least two atoms on the candidate ring have substituents that are parts
    of other benzene rings connected by a single bond (i.e. the candidate ring and the substituent ring 
    share exactly one atom), and whether these two substituents are in para positions (separated by three bonds along the cycle).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule has a para-terphenyl skeleton, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Error sanitizing molecule: {e}"
    
    ring_info = mol.GetRingInfo().AtomRings()  # list of tuples, each a set of atom indices in a ring

    # Helper: Check if the ring (given as tuple of atom indices) is a benzene ring.
    def is_benzene_ring(ring):
        if len(ring) != 6:
            return False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() != "C" or not atom.GetIsAromatic():
                return False
        return True

    # Get all benzene rings as sets.
    benzene_rings = [set(ring) for ring in ring_info if is_benzene_ring(ring)]
    if not benzene_rings:
        return False, "No benzene rings found"

    # Helper: Given a candidate ring (set of atom indices), determine its cyclic ordering.
    # Since every atom in a proper benzene ring has exactly two neighbors within the ring,
    # we can do a simple traversal starting from an arbitrary atom.
    def get_cycle_ordering(candidate):
        candidate = set(candidate)  # ensure it's a set
        # Pick an arbitrary starting atom.
        start = next(iter(candidate))
        ordering = [start]
        prev = None
        current = start
        # In a benzene ring there are 6 atoms.
        while len(ordering) < len(candidate):
            curr_atom = mol.GetAtomWithIdx(current)
            # Look for neighbors that are also in the candidate set and are not the previous atom.
            found = False
            for nbr in curr_atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx in candidate and nbr_idx != prev:
                    ordering.append(nbr_idx)
                    prev = current
                    current = nbr_idx
                    found = True
                    break
            if not found:
                # this should not happen for a proper ring
                break
        return ordering

    # Now, for each benzene ring candidate,
    # try to identify it as the central ring of a para-terphenyl skeleton.
    for candidate_ring in benzene_rings:
        # Reconstruct the cyclic ordering.
        ordering = get_cycle_ordering(candidate_ring)
        if len(ordering) != 6:
            continue  # not a proper benzene cycle
        
        # For each atom on the candidate ring, check for a substituent benzene ring.
        # A valid substituent is an aromatic benzene ring that shares exactly one atom with the candidate.
        substituent_positions = []
        for pos, atom_idx in enumerate(ordering):
            atom = mol.GetAtomWithIdx(atom_idx)
            valid_substituent_found = False
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                # Only consider if neighbor is NOT part of the candidate ring.
                if nbr_idx in candidate_ring:
                    continue
                # Look for any benzene ring that contains this neighbor.
                for benzene in benzene_rings:
                    # To be a substituent the benzene ring should share exactly one atom (the candidate atom)
                    if nbr_idx in benzene:
                        overlap = candidate_ring.intersection(benzene)
                        if len(overlap) == 1:
                            valid_substituent_found = True
                            break
                if valid_substituent_found:
                    substituent_positions.append(pos)
                    break  # Only one substituent per candidate atom
        # We need at least two substituents.
        if len(substituent_positions) < 2:
            continue  # try next candidate central ring
        
        # Check if any two substituents are para to each other.
        n = 6  # number of atoms in the ring
        for i in range(len(substituent_positions)):
            for j in range(i+1, len(substituent_positions)):
                diff = abs(substituent_positions[i] - substituent_positions[j])
                circular_diff = min(diff, n - diff)
                if circular_diff == 3:  # exactly para in benzene
                    return True, "Found para-terphenyl skeleton: central benzene ring with benzene substituents in para positions attached via single (non-fused) bonds."
    return False, "Could not find a central benzene ring with two para-substituted benzene rings (non-fused) required for para-terphenyl."

# Example usage:
# Uncomment below to test on sample SMILES:
# smiles_list = [
#     "O=C(OC1=C(OC(=O)C)C(=C(O)C(=C1C2=CC(O)=C(O)C=C2)O)C3=CC=C(O)C=C3)C",  # 2',3'-diacetoxy-... p-terphenyl
#     "COc1cc(-c2ccccc2)c(OC)c(O)c1-c1ccc(O)c(O)c1",  # 3-hydroxyterphenyllin, similar derivative
# ]
# for s in smiles_list:
#     result, reason = is_para_terphenyl(s)
#     print(f"SMILES: {s}\nResult: {result}, Reason: {reason}\n")