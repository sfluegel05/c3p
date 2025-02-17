"""
Classifies: CHEBI:75874 para-terphenyl
"""
"""
Classifies: para-terphenyl derivatives (based on the 1,4-diphenylbenzene skeleton).
Definition: A ring assembly based on a 1,4-diphenylbenzene skeleton and its substituted derivatives thereof.
"""

from rdkit import Chem

def is_para_terphenyl(smiles: str):
    """
    Determines if a molecule is a para-terphenyl derivative based on its SMILES string.
    A para-terphenyl has a central benzene ring with two benzene (phenyl) rings attached
    at positions that are opposite (para) on the central ring.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule has a para-terphenyl skeleton, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure aromaticity perception is done
    Chem.SanitizeMol(mol)
    
    ring_info = mol.GetRingInfo().AtomRings()  # list of tuples of atom indices in each ring
    
    # Helper function to check if a ring is a benzene (6-membered aromatic ring with all C)
    def is_benzene_ring(ring):
        if len(ring) != 6:
            return False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() != "C" or not atom.GetIsAromatic():
                return False
        return True
    
    # For each ring, try to see if it can be the central ring
    for ring in ring_info:
        if not is_benzene_ring(ring):
            continue  # skip rings that are not benzene
        
        # Attempt to get the ring atoms in a cyclic order.
        # The ring info from RDKit is in order (cyclic order).
        ring_atoms = list(ring)
        n = len(ring_atoms)  # Should be 6 for benzene
        
        # Create a map from atom index to its position in the ring order
        pos_in_ring = {atom_idx: i for i, atom_idx in enumerate(ring_atoms)}
        
        # For each atom in this central ring, check neighbors (outside the ring)
        substituent_positions = []
        for i, atom_idx in enumerate(ring_atoms):
            atom = mol.GetAtomWithIdx(atom_idx)
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx in ring:
                    continue  # ignore neighbors that are within the ring
                # Check if the neighbor atom belongs to a benzene ring (phenyl group)
                in_benzene = False
                # Check all rings that include the neighbor; if one is benzene then treat as phenyl substituent
                for ring_candidate in ring_info:
                    if nbr_idx in ring_candidate and is_benzene_ring(ring_candidate):
                        in_benzene = True
                        break
                if in_benzene:
                    substituent_positions.append(i)
                    # break after first found benzene substituent on that atom
                    break
        
        # For a para substitution, we need two substituents on the ring with positions separated by 3 atoms around the ring
        if len(substituent_positions) < 2:
            continue  # not enough substituents on this ring
        
        # Check every pair for circular distance of 3 (in a hexagon, para positions are 3 apart)
        for i in range(len(substituent_positions)):
            for j in range(i+1, len(substituent_positions)):
                diff = abs(substituent_positions[i] - substituent_positions[j])
                circular_diff = min(diff, n - diff)
                if circular_diff == 3:
                    return True, "Found para-terphenyl skeleton: central benzene ring with benzene substituents in para positions."
        
    return False, "Could not find a central benzene ring with two para-substituted benzene rings."

# Example usage (uncomment to test):
# smiles_examples = [
#     "O(C1=C(O)C(C2=CC=CC=C2)=CC(=C1C3=CC=CC=C3)O)C",  # 2-(2,5-dihydroxyphenyl)-3-methoxy-5-phenylbenzene-1,4-diol
#     "COc1cc(-c2ccccc2)c(OC)c(O)c1-c1ccc(O)cc1"          # 4''-deoxyterphenyllin
# ]
# for s in smiles_examples:
#     result, message = is_para_terphenyl(s)
#     print(f"SMILES: {s}\nResult: {result}, Reason: {message}\n")