"""
Classifies: CHEBI:37581 gamma-lactone
"""
#!/usr/bin/env python3
"""
Classifies: gamma-lactone
Definition: A gamma-lactone (for our purpose) is a cyclic ester that
  contains a five-membered ring made up of exactly four carbons and one oxygen,
  with one of the ring carbons (the lactone carbonyl carbon) bearing an exocyclic 
  carbonyl group (i.e. a double bond to an oxygen that is not part of the ring).
  
This code uses RDKit to parse the SMILES and then manually inspects each five-membered ring.
"""

from rdkit import Chem

def is_gamma_lactone(smiles: str):
    """
    Determines if a molecule is a gamma‐lactone based on its SMILES string.
    A gamma‐lactone (by our strict definition) must contain at least one five‐membered
    cyclic ester ring meeting these criteria:
      • The ring is composed of 5 atoms: exactly one oxygen and exactly four carbons.
      • One of the ring carbons (the lactone carbonyl carbon) has an exocyclic oxygen
        attached by a double bond. In addition, the exocyclic oxygen should be terminal 
        (only attached to that carbon), helping to ensure that it is a carbonyl.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): A tuple where the first element is True if a gamma‐lactone motif is 
                     found, and the second element gives an explanation; otherwise, False 
                     with a reason. If the SMILES is invalid, returns (False, "Invalid SMILES string").
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information (each ring is a tuple of atom indices)
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()

    # Iterate over each ring
    for ring in rings:
        # We only consider five-membered rings.
        if len(ring) != 5:
            continue
        
        # Count atoms in the ring by element and ensure that only carbon and oxygen appear.
        num_oxygens = 0
        num_carbons = 0
        valid_ring = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            Z = atom.GetAtomicNum()
            if Z == 8:
                num_oxygens += 1
            elif Z == 6:
                num_carbons += 1
            else:
                valid_ring = False  # any atom other than C or O disqualifies the ring under our strict definition
                break
        if not valid_ring:
            continue

        # We require exactly 1 oxygen and 4 carbons in the ring.
        if num_oxygens != 1 or num_carbons != 4:
            continue

        # Now check for the exocyclic carbonyl group.
        # For each carbon in the ring, check its neighbors that are not in the ring.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue  # only a carbon can be the carbonyl carbon
            ring_neighbors = {nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetIdx() in ring}
            exo_carbonyl_found = False
            for nbr in atom.GetNeighbors():
                # Only consider neighbors not in the ring.
                if nbr.GetIdx() in ring:
                    continue
                # Check that neighbor is oxygen and that the bond is a double bond.
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if nbr.GetAtomicNum() == 8 and bond.GetBondType() == Chem.BondType.DOUBLE:
                    # Additionally, check that this exocyclic oxygen is "terminal" (only attached to the carbon).
                    if nbr.GetDegree() == 1:
                        exo_carbonyl_found = True
                        break
            if exo_carbonyl_found:
                return True, ("Found a five‐membered lactone ring (4 carbons, 1 oxygen) with an exocyclic "
                              "carbonyl group attached to one of the ring carbons.")
    
    return False, "No gamma‐lactone ring (five‐membered cyclic ester with proper connectivity) found."

# Example usage (for testing purposes):
if __name__ == "__main__":
    # Test with a known gamma-lactone (gamma-nonanolactone)
    test_smiles = "CCCCCC[C@H]1CCC(=O)O1"
    result, reason = is_gamma_lactone(test_smiles)
    print("Result:", result)
    print("Reason:", reason)