"""
Classifies: CHEBI:37581 gamma-lactone
"""
#!/usr/bin/env python3
"""
Classifies: gamma-lactone
Definition: A lactone having a five‐membered lactone ring.
A gamma‐lactone (by our definition) is a cyclic ester formed by a five‐membered ring made up 
of exactly four carbons and one oxygen (the ring oxygen), where one of the ring carbons 
(called the lactone carbonyl) carries an exocyclic carbonyl (C=O) group (with the oxygen being 
external to the ring). This approach explicitly examines each five‐membered ring in the molecule.
"""

from rdkit import Chem

def is_gamma_lactone(smiles: str):
    """
    Determines if a molecule is a gamma‐lactone based on its SMILES string.
    A gamma‐lactone (as defined here) is a molecule that contains at least one five‐membered cyclic ester
    which meets the following criteria:
      • The ring is composed of exactly five atoms.
      • Exactly one of the ring atoms is an oxygen (the ring oxygen) and the other four are carbons.
      • One of those carbons (the lactone carbonyl carbon) is attached (via a double bond)
        to an oxygen that is not part of the ring.
    While many molecules may contain a five‐membered lactone ring, stereochemistry and additional fused rings 
    may complicate the recognition. This method tries to ensure that the match is consistent with a gamma‐lactone.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): A tuple where the first element is True if a gamma‐lactone motif is found, 
                     and the second element gives an explanation; otherwise False with a reason.
                     If the SMILES is invalid, returns (False, "Invalid SMILES string").
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information from the molecule
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # each ring is a tuple of atom indices

    for ring in rings:
        if len(ring) != 5:
            continue  # Only interested in 5-membered rings

        # Count ring atoms by element:
        num_oxygens = 0
        num_carbons = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            Z = atom.GetAtomicNum()
            if Z == 8:
                num_oxygens +=1
            elif Z == 6:
                num_carbons +=1
            else:
                # any other element disqualifies this ring per our strict definition
                num_carbons = -999  # force failure
                break
        # Our definition requires exactly one oxygen and exactly four carbons.
        if num_oxygens != 1 or num_carbons != 4:
            continue

        # Now check if one of the carbons in this ring has an exocyclic carbonyl group.
        # That means: it must have a double bond to an oxygen that is not part of the ring.
        lactone_carbon_found = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                # Only carbons can be the carbonyl carbon
                continue
            # Check neighbors of this carbon:
            for neighbor, bond in [(nbr, mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx()))
                                     for nbr in atom.GetNeighbors()]:
                # Skip neighbors that are part of the ring.
                if neighbor.GetIdx() in ring:
                    continue
                # Look for a double bond to an oxygen
                if neighbor.GetAtomicNum() == 8 and bond.GetBondType() == Chem.BondType.DOUBLE:
                    lactone_carbon_found = True
                    break
            if lactone_carbon_found:
                break

        if lactone_carbon_found:
            return True, ("Found a five‐membered lactone ring composed of four carbons and one oxygen "
                          "with an exocyclic carbonyl group attached to one of the ring carbons.")
    
    # If we exit the loop, no gamma‐lactone ring has been identified.
    return False, "No gamma‐lactone ring (five‐membered cyclic ester with proper connectivity) found."

# Example usage (for testing purposes):
if __name__ == "__main__":
    # Example: gamma-nonanolactone SMILES from the provided list.
    test_smiles = "CCCCCC[C@H]1CCC(=O)O1"
    result, reason = is_gamma_lactone(test_smiles)
    print("Result:", result)
    print("Reason:", reason)