"""
Classifies: CHEBI:50753 isoflavonoid
"""
"""
Classifies: Any 1-benzopyran with an aryl substituent at position 3 (isoflavonoid).
An isoflavonoid is defined as a benzopyran (a 6‐membered pyran ring containing exactly one oxygen)
fused to a benzene ring (6‐membered, fully aromatic, carbon only) such that an exocyclic aryl substituent
(i.e. a 6‐membered aromatic ring not part of the fused core) is attached at the “3‐position” of the pyran.
This code implements a heuristic that uses ring order (via RDKit’s GetSymmSSSR) so that the position
two bonds away (clockwise from the oxygen in the ordered pyran) is taken as position 3.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid based on its SMILES string.
    An isoflavonoid is defined as any 1-benzopyran with an aryl substituent at position 3.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an isoflavonoid, otherwise False.
        str: Explanation of the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"

    # Get basic ring information from the molecule.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # each ring is a tuple of atom indices
    if not atom_rings:
        return False, "No rings found in the molecule"

    # Use RDKit's Symmetry-based ring finder for ordered rings.
    # This returns RingInfo objects (as tuples representing a cyclic order)
    rings = list(Chem.GetSymmSSSR(mol))
    
    # Build candidate benzene rings: 6-membered rings with all aromatic carbons.
    benzene_rings = []
    for ring in atom_rings:
        if len(ring) == 6:
            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() and
                   mol.GetAtomWithIdx(idx).GetAtomicNum() == 6
                   for idx in ring):
                benzene_rings.append(set(ring))
    if not benzene_rings:
        return False, "No benzene rings (6-membered, fully aromatic carbons) found."

    # Collect candidate pyran rings from the ordered rings (from GetSymmSSSR):
    # We require a 6-membered ring with exactly one oxygen and five carbons.
    pyran_candidates = []  # list of tuples (ordered_ring, set(ring))
    for cyc in rings:
        if len(cyc) == 6:
            atoms = [mol.GetAtomWithIdx(idx) for idx in cyc]
            ocount = sum(1 for a in atoms if a.GetAtomicNum() == 8)
            ccount = sum(1 for a in atoms if a.GetAtomicNum() == 6)
            if ocount == 1 and ccount == 5:
                pyran_candidates.append( (list(cyc), set(cyc)) )
    if not pyran_candidates:
        return False, "No candidate pyran ring (6-membered with exactly one oxygen) found."
    
    # Now try to detect a fused benzopyran system with the correct exocyclic substitution.
    for ordered_pyran, pyran_set in pyran_candidates:
        # Look for a fused benzene ring which shares 2 or 3 atoms with the pyran.
        for benzene in benzene_rings:
            shared = pyran_set.intersection(benzene)
            if len(shared) in (2, 3):
                fused_core = pyran_set.union(benzene)
                # Try to locate the oxygen atom in the ordered pyran ring.
                oxygen_index = None
                for i, idx in enumerate(ordered_pyran):
                    atom = mol.GetAtomWithIdx(idx)
                    if atom.GetAtomicNum() == 8:
                        oxygen_index = i
                        break
                if oxygen_index is None:
                    continue  # should not occur as we filtered for one oxygen

                # In a typical 1-benzopyran, the heterocyclic numbering (starting at the oxygen as 1)
                # implies that the atom two positions away (clockwise) is at position 3.
                candidate_pos = ordered_pyran[(oxygen_index + 2) % 6]
                # If this candidate atom is already in the fused interface, skip it.
                if candidate_pos in shared:
                    # It would be fused with the benzene, not free for exocyclic substitution.
                    continue

                candidate_atom = mol.GetAtomWithIdx(candidate_pos)
                # Look for neighbors of candidate_atom that are not in the pyran ring.
                for nbr in candidate_atom.GetNeighbors():
                    if nbr.GetIdx() in pyran_set:
                        continue  # Skip connections within the pyran
                    # For an aryl substituent, the neighbor should be an aromatic carbon.
                    if nbr.GetAtomicNum() != 6 or not nbr.GetIsAromatic():
                        continue
                    # Now check that this neighbor participates in a 6-membered aromatic ring (phenyl-like)
                    found_aryl_ring = False
                    for ring in atom_rings:
                        if len(ring) == 6 and nbr.GetIdx() in ring:
                            if all(mol.GetAtomWithIdx(i).GetIsAromatic() and mol.GetAtomWithIdx(i).GetAtomicNum() == 6 
                                   for i in ring):
                                # Ensure this ring is not completely part of the fused benzopyran core.
                                if not set(ring).issubset(fused_core):
                                    found_aryl_ring = True
                                    break
                    if found_aryl_ring:
                        return True, ("Molecule contains a fused benzopyran core (6-membered pyran with one oxygen fused "
                                      "with a benzene ring sharing 2–3 atoms) and an exocyclic aryl substituent attached at "
                                      "the position expected to be position 3.")
    # If we have exhausted the candidates without a positive identification.
    return False, "Scaffold not recognized as isoflavonoid (no fused benzopyran core with an aryl substituent at position 3 found)."

# Example usage (for testing purposes – you may run several SMILES strings):
if __name__ == '__main__':
    test_smiles = [
        # Some test SMILES examples (these do not cover full complexity)
        "O1C(C2=CC=CC=C2)=CC(O)=C1",  # simple benzopyran (may lack the proper substitution)
        "COc1ccc(-c2coc3ccccc3c2)cc1"  # a rough example of benzopyran with a phenyl substituent.
    ]
    for smi in test_smiles:
        res, reason = is_isoflavonoid(smi)
        print(f"SMILES: {smi}\nResult: {res}\nReason: {reason}\n")