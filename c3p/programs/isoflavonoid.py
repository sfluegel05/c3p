"""
Classifies: CHEBI:50753 isoflavonoid
"""
"""
Classifies: Any 1-benzopyran with an aryl substituent at position 3 (isoflavonoid).
The classifier uses a heuristic that attempts to identify a pyran ring (six‐membered ring with one oxygen)
that is fused (shares ≥2 atoms) with a benzene ring. Then it looks for an exocyclic aromatic (aryl) substituent
attached to one of the non-fused atoms of the pyran ring.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid based on its SMILES string.
    An isoflavonoid is defined as any 1-benzopyran with an aryl substituent at position 3.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as an isoflavonoid, otherwise False.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring info (each ring as tuple of atom indices)
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found in the molecule"

    # Collect candidate benzene rings (aromatic, size 6)
    benzene_rings = []
    for ring in ring_info:
        if len(ring) == 6:
            # Check if all atoms in the ring are aromatic and are carbons
            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() and mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
                benzene_rings.append(set(ring))
    if not benzene_rings:
        return False, "No benzene rings (aromatic 6-membered all-carbon) found for a fused system"

    # Now, look for candidate pyran rings:
    # These are 6-membered rings that contain exactly one oxygen.
    candidate_pyran = []
    for ring in ring_info:
        if len(ring) == 6:
            atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            oxy_count = sum(1 for a in atoms if a.GetAtomicNum() == 8)
            # Candidate pyran: exactly one oxygen.
            if oxy_count == 1:
                candidate_pyran.append(set(ring))
    if not candidate_pyran:
        return False, "No six-membered ring with one oxygen (pyran-like ring) found"

    # For each candidate pyran, check if it is fused with one of the benzene rings.
    for pyran in candidate_pyran:
        fusion_found = False
        shared_ring = None
        for benzene in benzene_rings:
            # Fusion typically means sharing at least 2 atoms.
            if len(pyran.intersection(benzene)) >= 2:
                fusion_found = True
                shared_ring = benzene
                break
        if not fusion_found:
            continue  # try next candidate pyran

        # In a 1-benzopyran, the pyran ring (sometimes called the C-ring of flavonoids)
        # is fused with a benzene ring (the A-ring). The aryl substituent (B-ring) is attached
        # at the position that is not participating in the fusion (typically position 3).
        non_fused_atoms = pyran.difference(shared_ring)
        # Heuristically, look for an atom in the candidate pyran (outside the fusion)
        # that has a neighbor not in the pyran and that neighbor is part of an aromatic ring.
        for atom_idx in non_fused_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            # Consider neighbors that are not in the pyran core (i.e. exocyclic)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in pyran:
                    continue
                # Check if the neighbor is aromatic and is carbon-based
                if nbr.GetAtomicNum() != 6 or not nbr.GetIsAromatic():
                    continue
                # Now check if this neighbor is part of a 6-membered aromatic ring (phenyl group)
                nbr_in_ring = False
                for ring in ring_info:
                    if len(ring) == 6 and nbr.GetIdx() in ring:
                        # Verify that all atoms in the ring are aromatic carbons
                        if all(mol.GetAtomWithIdx(i).GetAtomicNum() == 6 and mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
                            nbr_in_ring = True
                            break
                if nbr_in_ring:
                    # We found an aryl substituent attached to an atom of the candidate pyran.
                    # (It is difficult to unambiguously assign the number '3' without full standard numbering,
                    #  but the heuristics strongly suggest the desired motif.)
                    return True, "Molecule contains a benzopyran core (fused pyran and benzene rings) with an exocyclic aryl substituent."
    return False, "Scaffold not recognized as isoflavonoid (no benzopyran core with aryl substituent found)."

# Example usage:
if __name__ == '__main__':
    # A couple of simple tests (users can try with the provided SMILES examples for isoflavonoids)
    test_smiles = [
        "O1C(C2=CC=CC=C2)=CC(O)=C1",  # a simple benzopyran – may or may not have an aryl substituent
        "COc1ccc(-c2coc3ccccc3c2)cc1"  # a rough example showing an aryl group attached to an oxygen heterocycle
    ]
    for smi in test_smiles:
        result, reason = is_isoflavonoid(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")