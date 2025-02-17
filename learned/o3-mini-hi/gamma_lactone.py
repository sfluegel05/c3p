"""
Classifies: CHEBI:37581 gamma-lactone
"""
#!/usr/bin/env python3
"""
Classifies: gamma-lactone
Definition: A lactone having a five-membered lactone ring.
In a gamma‐lactone the carbonyl carbon is part of a five‐membered ring that contains
exactly one oxygen atom (the ester oxygen) and the carbonyl group (C=O) has its oxygen
external to the ring.
"""

from rdkit import Chem

def is_gamma_lactone(smiles: str):
    """
    Determines if a molecule is a gamma-lactone based on its SMILES string.
    A gamma-lactone is a cyclic ester formed by a five-membered ring that includes exactly
    one ring oxygen. In this ring, the carbonyl carbon (a ring carbon) must have a double bond
    to an oxygen that is not part of the ring (the carbonyl oxygen) and be directly connected 
    to the ring oxygen (forming the ester linkage).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule qualifies as a gamma-lactone, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get all defined rings in the molecule.
    ring_atom_indices = mol.GetRingInfo().AtomRings()
    
    # Loop over rings and only consider those with exactly five atoms.
    for ring in ring_atom_indices:
        if len(ring) != 5:
            continue

        # Count how many oxygen atoms are in this ring.
        oxygen_in_ring = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        # In a typical gamma-lactone the ring itself (the cyclic ester) has exactly one oxygen.
        if oxygen_in_ring != 1:
            continue

        # Look for a carbon in the ring that shows the lactone connectivity.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Consider only carbons as the potential carbonyl carbon.
            if atom.GetAtomicNum() != 6:
                continue
            
            # Check that this carbon is directly bonded to the ring oxygen.
            has_ring_oxygen_neighbor = any(neighbor.GetIdx() in ring and neighbor.GetAtomicNum() == 8 
                                           for neighbor in atom.GetNeighbors())
            if not has_ring_oxygen_neighbor:
                continue

            # Now check that the same carbon has a double bond to an oxygen that is NOT in the ring.
            has_exterior_carbonyl = False
            for bond in atom.GetBonds():
                # Look for a double bond
                if bond.GetBondType() != Chem.BondType.DOUBLE:
                    continue
                other = bond.GetOtherAtom(atom)
                # We want a double-bonded oxygen that is outside of the ring.
                if other.GetAtomicNum() == 8 and other.GetIdx() not in ring:
                    has_exterior_carbonyl = True
                    break
            
            if has_exterior_carbonyl:
                return True, "Found five-membered lactone ring with ester connectivity (carbonyl outside, ring oxygen inside)."
    
    # If we cannot find any ring satisfying our criteria, return False.
    return False, "No gamma-lactone ring (five-membered cyclic ester with proper connectivity) found."

# Example usage (for testing purposes):
if __name__ == "__main__":
    # Test an example gamma-lactone: gamma-nonanolactone
    test_smiles = "CCCCCC[C@H]1CCC(=O)O1"
    result, reason = is_gamma_lactone(test_smiles)
    print("Result:", result)
    print("Reason:", reason)