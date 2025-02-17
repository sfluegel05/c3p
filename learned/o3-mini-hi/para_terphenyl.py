"""
Classifies: CHEBI:75874 para-terphenyl
"""
"""
Classifies: para-terphenyl derivatives
Definition: A ring assembly based on a 1,4-diphenylbenzene skeleton and its substituted derivatives thereof.
A true para-terphenyl contains a central benzene ring (non-fused) that is para-substituted 
(i.e. substituents at positions 1 and 4; in a hexagon these are 3 atoms apart) 
with benzene rings attached via a single (non-fused) bond.
"""

from rdkit import Chem

def is_para_terphenyl(smiles: str):
    """
    Determines if a molecule is a para-terphenyl derivative based on its SMILES string.
    
    The algorithm finds a benzene ring (6-membered aromatic carbocycle) 
    in which at least two substituents (attached by a single bond) are themselves 
    benzene rings attached with exactly one overlapping atom with the central ring.
    These substituents should lie para relative to each other (i.e. positions 0 and 3 in the ring order).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule has a para-terphenyl skeleton, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure aromaticity and ring info are properly set.
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Error sanitizing molecule: {e}"
    
    ring_info = mol.GetRingInfo().AtomRings()  # list of tuples; each tuple is indices of atoms in a ring

    # Helper: Check if a ring (atom indices) is a benzene ring.
    def is_benzene_ring(ring):
        if len(ring) != 6:
            return False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() != "C" or not atom.GetIsAromatic():
                return False
        return True

    # List all benzene rings in the molecule.
    benzene_rings = [set(ring) for ring in ring_info if is_benzene_ring(ring)]
    if not benzene_rings:
        return False, "No benzene rings found"

    # For each benzene ring, test if it can be the central ring of a para-terphenyl.
    for candidate_ring in benzene_rings:
        candidate_list = list(candidate_ring)
        if len(candidate_list) != 6:
            continue  # not possible, but safe-check

        # Map each atom in candidate_ring to its position in the candidate list.
        # The ring info from RDKit is in cyclic order.
        pos_in_ring = {atom_idx: i for i, atom_idx in enumerate(candidate_list)}
        
        # Record positions on the central ring that have valid benzene substituents.
        substituent_positions = []
        
        # For each atom in the candidate central ring...
        for pos, atom_idx in enumerate(candidate_list):
            atom = mol.GetAtomWithIdx(atom_idx)
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                # Only consider substituents not part of the central ring.
                if nbr_idx in candidate_ring:
                    continue
                
                # Look for a benzene ring that the neighbor belongs to.
                # Also check that the overlap between that benzene ring and the candidate ring is exactly one.
                valid_substituent = False
                for benzene in benzene_rings:
                    if nbr_idx in benzene:
                        # Calculate the overlap with the candidate central ring.
                        overlap = candidate_ring.intersection(benzene)
                        if len(overlap) == 1:
                            valid_substituent = True
                            break
                if valid_substituent:
                    substituent_positions.append(pos)
                    # Only count one substituent per atom on the candidate ring.
                    break
        
        # If there are at least two substituents, check if any pair are para (3 bonds apart in a 6-membered ring)
        if len(substituent_positions) < 2:
            continue  # try next candidate central ring
        
        n = 6  # for benzene ring
        # Check every unique pair of substituent positions for para relationship.
        for i in range(len(substituent_positions)):
            for j in range(i+1, len(substituent_positions)):
                diff = abs(substituent_positions[i] - substituent_positions[j])
                circular_diff = min(diff, n - diff)
                if circular_diff == 3:
                    return True, "Found para-terphenyl skeleton: central benzene ring with benzene substituents in para positions attached via single (non-fused) bonds."
    return False, "Could not find a central benzene ring with two para-substituted benzene rings (non-fused) required for para-terphenyl."
    
# Example usage (uncomment to test):
# smiles_list = [
#     "O=C(OC1=C(OC(=O)C)C(=C(O)C(=C1C2=CC(O)=C(O)C=C2)O)C3=CC=C(O)C=C3)C",  # 2',3'-diacetoxy-... p-terphenyl, true positive
#     "CC(=O)C1=CC2=C(C=C1)C3=CC=CC=C3C=C2",  # 1-(2-phenanthrenyl)ethanone, false positive in previous algorithm
# ]
# for s in smiles_list:
#     result, reason = is_para_terphenyl(s)
#     print(f"SMILES: {s}\nResult: {result}, Reason: {reason}\n")