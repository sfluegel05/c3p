"""
Classifies: CHEBI:140310 phenyl acetates
"""
"""
Classifies: Phenyl Acetates
Definition: An acetate ester obtained by formal condensation of the carboxy group of acetic acid
with the hydroxy group of any phenol. For our purposes the molecule must contain an acetoxy group –O–C(=O)CH3
that is acyclic (i.e. the O–C(=O) bond is not part of a ring) and where the oxygen (O) is attached to an aromatic carbon
that is part of a benzene ring (a six-membered ring composed solely of aromatic carbons).
"""

from rdkit import Chem

def is_phenyl_acetates(smiles: str):
    """
    Determines whether a molecule qualifies as a phenyl acetate based on its SMILES string.

    The molecule must contain an acyclic acetoxy group (–O–C(=O)CH3) such that the oxygen (O) attached to acyl group
    is directly connected to an aromatic carbon that is a member of a benzene ring.

    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as a phenyl acetate, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS for an acetoxy group: [O][C](=O)[CH3]
    # We do not require the attached O to be aromatic here.
    acetoxy_smarts = "[O:1][C:2](=O)[CH3]"
    acetoxy_pattern = Chem.MolFromSmarts(acetoxy_smarts)
    if acetoxy_pattern is None:
        return False, "Error creating SMARTS pattern for acetoxy group"
    
    matches = mol.GetSubstructMatches(acetoxy_pattern)
    if not matches:
        return False, "No acetoxy group (-O-C(=O)CH3) found in the molecule"
    
    # Get ring information once for efficiency.
    ring_info = mol.GetRingInfo().AtomRings()
    
    # Loop over each acetoxy match candidate.
    for match in matches:
        # match is a tuple of atom indices corresponding to: O, carbonyl C, methyl C
        if len(match) != 3:
            continue  # we expect exactly three atoms
        
        o_idx, carbonyl_idx, methyl_idx = match
        
        o_atom = mol.GetAtomWithIdx(o_idx)
        # Ensure the acetoxy group is acyclic:
        # Check if the O--C(carbonyl) bond is in any ring.
        in_ring = False
        for ring in ring_info:
            if o_idx in ring and carbonyl_idx in ring:
                in_ring = True
                break
        if in_ring:
            # Skip this match because it is likely part of a lactone.
            continue

        # Find the neighbor of the oxygen atom that is not the carbonyl carbon.
        # In an ester, oxygen should have two neighbors: one carbonyl carbon and the other is the substituent.
        o_neighbors = [nbr for nbr in o_atom.GetNeighbors() if nbr.GetIdx() != carbonyl_idx]
        if not o_neighbors:
            continue  # no substituent attached?
        
        # Check if any neighbor is aromatic and part of a benzene ring.
        aromatic_aromatic_connection = False
        for nbr in o_neighbors:
            if not nbr.GetIsAromatic():
                continue  # not aromatic, skip.
            nbr_idx = nbr.GetIdx()
            # Look for a benzene ring: exactly 6 atoms, all aromatic carbons.
            for ring in ring_info:
                if nbr_idx in ring and len(ring) == 6:
                    # Check that each atom in the ring is a carbon and is aromatic.
                    if all(mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 6 and 
                           mol.GetAtomWithIdx(atom_idx).GetIsAromatic() for atom_idx in ring):
                        aromatic_aromatic_connection = True
                        break
            if aromatic_aromatic_connection:
                break
                
        if aromatic_aromatic_connection:
            return True, "Molecule contains an acyclic acetoxy group (-O-C(=O)CH3) attached to a benzene ring."
    
    # If none of the acetoxy groups qualify, provide a reason.
    return False, "Acetoxy group found but did not meet phenyl acetate criteria (either in a ring or not attached to a benzene ring)."

# Example usage:
if __name__ == "__main__":
    # A known phenyl acetate example: phenyl acetate.
    test_smiles = "CC(=O)Oc1ccccc1"  
    result, reason = is_phenyl_acetates(test_smiles)
    print("Test (phenyl acetate):", result, reason)
    
    # Testing one of the provided examples:
    test_smiles_example = "COc1cc(ccc1OC(C)=O)C(C)=O"  # 4-acetyl-2-methoxyphenyl acetate
    result, reason = is_phenyl_acetates(test_smiles_example)
    print("Test (4-acetyl-2-methoxyphenyl acetate):", result, reason)