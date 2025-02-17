"""
Classifies: CHEBI:38958 indole alkaloid
"""
#!/usr/bin/env python
"""
Classifies: Indole Alkaloid - An alkaloid containing an indole skeleton.
Improved approach: Instead of relying solely on a single SMARTS pattern, we
examine the ring system. In an indole, an aromatic five-membered ring (with one nitrogen)
is fused with an aromatic six-membered ring (sharing exactly two atoms). This procedure
helps capture variations of the indole core present in many alkaloids.
"""

from rdkit import Chem

def is_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is an indole alkaloid based on its SMILES string.
    An indole alkaloid is defined as an alkaloid that contains an indole skeleton,
    namely a bicyclic system consisting of a benzene ring fused with a pyrrole ring.
    
    This function first parses the SMILES string. It then examines the ring systems:
    it looks for at least one aromatic five-membered ring that contains exactly one nitrogen atom
    and a fused aromatic six-membered ring (sharing exactly two atoms with the 5-membered ring).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule contains an indole skeleton, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string to an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure aromaticity is computed.
    Chem.SanitizeMol(mol)
    
    # Get all rings in the molecule as tuples of atom indices.
    ring_info = mol.GetRingInfo()
    if not ring_info:
        return False, "No rings found in the molecule"
    rings = ring_info.AtomRings()
    
    # Helper function: check if all atoms in a ring are aromatic.
    def is_aromatic_ring(atom_indices):
        return all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in atom_indices)
    
    # Separate rings by size and type.
    fm_rings = []  # five-membered rings with one nitrogen (potential pyrrole rings)
    sm_rings = []  # six-membered rings (potential benzene rings)
    
    for ring in rings:
        if len(ring) == 5 and is_aromatic_ring(ring):
            # Count nitrogen atoms in the ring.
            n_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
            # In an indole, the five-membered ring should have exactly one nitrogen.
            if n_count == 1:
                fm_rings.append(set(ring))
        elif len(ring) == 6 and is_aromatic_ring(ring):
            sm_rings.append(set(ring))
    
    # Look for a pair (five-membered ring, six-membered ring) that share exactly two atoms.
    for fm in fm_rings:
        for sm in sm_rings:
            if len(fm.intersection(sm)) == 2:
                return True, "Molecule contains a fused aromatic five-membered ring (with one nitrogen) and six-membered ring characteristic of an indole skeleton."
    
    # If no matching fused ring pair is found, classification fails.
    return False, "Molecule does not contain the fused ring system characteristic of an indole skeleton."

# Example usage:
# test_smiles = "CN[C@@H]1C[C@H]2O[C@@](C)([C@@H]1OC)N1C3=C(C=CC=C3)C3=C1C1=C(C4=C(C=CC=C4)N21)C1=C3CNC1=O"  # staurosporine
# result, reason = is_indole_alkaloid(test_smiles)
# print(result, reason)