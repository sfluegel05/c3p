"""
Classifies: CHEBI:37581 gamma-lactone
"""
#!/usr/bin/env python3
"""
Classifies: gamma-lactone
Definition: For our purposes a gamma-lactone is a cyclic ester that contains a five-membered ring 
  made up of exactly four carbons and one oxygen, where one of the ring carbons 
  (the lactone carbonyl carbon) bears an exocyclic terminal carbonyl group (a C=O where the O is not in the ring).
  
This implementation uses RDKit’s symmetric SSSR finder to get all rings (even those in fused systems)
and then inspects each five-membered ring for the proper composition and connectivity.
"""

from rdkit import Chem

def is_gamma_lactone(smiles: str):
    """
    Determines if a molecule is a gamma‐lactone based on its SMILES string by looking for a five‐membered cyclic ester 
    ring composed of 4 carbons and 1 oxygen in which one of the carbons has a terminal exocyclic carbonyl group.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
       (bool, str): Tuple where the first element is True if a gamma-lactone motif is found, and the second element 
                    provides an explanation. Returns (False, reason) if not found, or if the SMILES string is invalid.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Use the symmetric SSSR finder to obtain all rings (even those in fused systems)
    rings = list(Chem.GetSymmSSSR(mol))
    
    # Iterate over each ring
    for ring in rings:
        ring_indices = set(ring)  # ring is a tuple of atom indices
        if len(ring_indices) != 5:
            continue  # only consider five-membered rings

        # Count number of carbon and oxygen atoms in the ring
        num_carbons = 0
        num_oxygens = 0
        for idx in ring_indices:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6:
                num_carbons += 1
            elif atom.GetAtomicNum() == 8:
                num_oxygens += 1
            else:
                # Ring has an atom other than carbon or oxygen–not our target ring.
                num_carbons = -1
                break
        if num_carbons == -1 or num_oxygens != 1 or num_carbons != 4:
            continue  # does not match gamma-lactone ring composition

        # Now check for the exocyclic carbonyl feature.
        # For each carbon atom in the ring, we check two things:
        # 1. It must be connected (within the ring) to the ring oxygen.
        # 2. It must have a neighbor (not in the ring) that is an oxygen connected by a double bond,
        #    and that oxygen must be terminal (degree 1).
        for idx in ring_indices:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue  # only carbon atoms can function as the lactone carbonyl carbon.
            # Check if one of its in-ring neighbors is the ring oxygen.
            in_ring_oxygens = [nbr for nbr in atom.GetNeighbors() if nbr.GetIdx() in ring_indices and nbr.GetAtomicNum() == 8]
            if not in_ring_oxygens:
                continue  # if no oxygen neighbor in the ring, skip
            # Now check for an exocyclic double-bonded terminal oxygen.
            exo_carbonyl_found = False
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring_indices:
                    continue  # skip atoms that are part of the ring
                if nbr.GetAtomicNum() == 8:
                    # Get bond between atom and neighbor
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                    if bond is None:
                        continue
                    # Check if bond order is double; note: sometimes aromatic bonds need special handling,
                    # but for carbonyl we expect an explicit double bond.
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        # Check that the exocyclic oxygen is terminal
                        if nbr.GetDegree() == 1:
                            exo_carbonyl_found = True
                            break
            if exo_carbonyl_found:
                return True, ("Found a five‐membered lactone ring (4 carbons, 1 oxygen) where a ring carbon is bonded to a terminal "
                              "exocyclic carbonyl oxygen.")
    
    return False, "No gamma‐lactone ring (five‐membered cyclic ester with proper connectivity) found."

# Example usage:
if __name__ == "__main__":
    # Test with a known gamma-lactone (gamma-nonanolactone)
    test_smiles = "CCCCCC[C@H]1CCC(=O)O1"
    result, reason = is_gamma_lactone(test_smiles)
    print("Result:", result)
    print("Reason:", reason)