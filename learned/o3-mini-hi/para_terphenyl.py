"""
Classifies: CHEBI:75874 para-terphenyl
"""
"""
Classifies: para-terphenyl derivatives
Definition: A ring assembly based on a 1,4-diphenylbenzene skeleton and its substituted derivatives.
A true para-terphenyl should possess a central benzene ring (six-membered, aromatic, composed of carbons)
with two benzene rings attached via a single (non-fused) bond at positions that are para (separated by 3 atoms).
"""

from rdkit import Chem

def is_para_terphenyl(smiles: str):
    """
    Determines if the molecule defined by a SMILES string exhibits a para-terphenyl motif.
  
    The method is as follows:
    1. Parse and sanitize the molecule.
    2. Retrieve all rings using the ring information from RDKit.
    3. Identify all six-membered rings in which every atom is aromatic carbon.
       These are treated both as candidates for central rings and external benzene substituents.
    4. For each candidate central ring, generate an ordered cycle (by traversing its connectivity).
    5. For each atom in the ordered ring, check its neighbors (outside of the ring) to see if any
       one belongs to an external benzene ring (non-fused, i.e. the intersection of the rings is exactly one atom).
    6. If two such substituents are found and the difference between their positions on the ring is 3
       (modulo 6), we classify the molecule as a para-terphenyl derivative.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if a para-terphenyl motif is found, else False.
        str: Explanation of the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Error sanitizing molecule: {e}"
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # tuple of tuples of atom indices
    
    # Helper: check if a ring (list of atom indices) is an aromatic benzene ring:
    def is_benzene_ring(ring):
        if len(ring) != 6:
            return False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if not (atom.GetIsAromatic() and atom.GetSymbol() == "C"):
                return False
        return True
    
    # Get all six-membered aromatic carbon rings
    benzene_rings = [list(r) for r in atom_rings if is_benzene_ring(r)]
    
    if not benzene_rings:
        return False, "No aromatic six-membered (benzene) rings found in the molecule."
    
    # Function to order the atoms in a ring by connectivity.
    # Assumes the ring is a simple cycle (as in benzene).
    def order_ring(ring):
        # ring is a list of atom indices; we build an ordered list by following connectivity.
        ordered = [ring[0]]
        # get neighbors within the ring for the starting atom
        current = ring[0]
        # find one neighbor that is in the ring
        neighs = [nbr.GetIdx() for nbr in mol.GetAtomWithIdx(current).GetNeighbors() if nbr.GetIdx() in ring]
        if not neighs:
            return ring  # fallback
        prev = current
        current = neighs[0]
        ordered.append(current)
        while len(ordered) < len(ring):
            atom_obj = mol.GetAtomWithIdx(current)
            # from neighbors in ring, pick one that is not the previous atom
            next_candidates = [nbr.GetIdx() for nbr in atom_obj.GetNeighbors() if (nbr.GetIdx() in ring and nbr.GetIdx() != prev)]
            if not next_candidates:
                break
            prev, current = current, next_candidates[0]
            ordered.append(current)
        return ordered

    # Now we want to try each benzene ring as a candidate central ring
    for central_ring in benzene_rings:
        ordered_central = order_ring(central_ring)
        substituent_positions = []  # positions on the central ring (0-indexed in the ordered list) where substituents attach
        
        n = 6  # number of atoms in a benzene ring (central)
        # For each atom in the ordered central ring, examine neighbors not in the ring.
        for pos, atom_idx in enumerate(ordered_central):
            central_atom = mol.GetAtomWithIdx(atom_idx)
            found_substituent = False
            for nbr in central_atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx in central_ring:
                    continue  # ignore atoms within the central ring
                # Check if this neighbor belongs to any benzene ring (other than the central ring) in a non-fused manner.
                for ext_ring in benzene_rings:
                    # Skip if the external ring is exactly the same as the central ring.
                    if set(ext_ring) == set(central_ring):
                        continue
                    if nbr_idx in ext_ring:
                        # A non-fused connection means that the intersection of the rings is exactly one atom.
                        if len(set(central_ring).intersection(ext_ring)) == 1:
                            substituent_positions.append(pos)
                            found_substituent = True
                            break
                if found_substituent:
                    break
        
        if len(substituent_positions) < 2:
            continue  # try next candidate central ring
        
        # Check whether any pair of substituents are in para positions (separated by 3 atoms in the 6-membered cycle).
        for i in range(len(substituent_positions)):
            for j in range(i+1, len(substituent_positions)):
                diff = abs(substituent_positions[i] - substituent_positions[j])
                circular_diff = min(diff, n - diff)
                if circular_diff == 3:
                    return True, ("Found para-terphenyl motif: a central benzene ring with two non-fused benzene substituents "
                                  f"attached at positions {substituent_positions[i]} and {substituent_positions[j]} (para to each other).")
    
    return False, "Could not find a central benzene ring with two para-attached benzene substituents (non-fused)."

# Example usage:
# Uncomment below to test with some SMILES strings:
# test_smiles = [
#    "O=C(OC1=C(OC(=O)C)C(=C(O)C(=C1C2=CC(O)=C(O)C=C2)O)C3=CC=C(O)C=C3)C",  # 2',3'-diacetoxy-3,4,5',6',4''-pentahydroxy-p-terphenyl 
#    "COc1cc(-c2ccccc2)c(OC)c(O)c1-c1ccc(O)c(O)c1",  # A simplified para-terphenyl derivative
# ]
# for s in test_smiles:
#     res, msg = is_para_terphenyl(s)
#     print(f"SMILES: {s}\nResult: {res}, Reason: {msg}\n")