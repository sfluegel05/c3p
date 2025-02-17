"""
Classifies: CHEBI:37581 gamma-lactone
"""
#!/usr/bin/env python3
"""
Classifies: gamma-lactone
Definition: A lactone having a five-membered lactone ring.
"""

from rdkit import Chem

def is_gamma_lactone(smiles: str):
    """
    Determines if a molecule is a gamma-lactone based on its SMILES string.
    A gamma-lactone is defined as a lactone that contains a five-membered ring in which 
    one carbon (the carbonyl carbon) is double-bonded to an oxygen (outside the ring) and 
    is also bonded to an oxygen that is a member of the ring (the ester oxygen).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule qualifies as a gamma-lactone, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Retrieve information about all rings in the molecule
    ring_info = mol.GetRingInfo().AtomRings()
    
    # Iterate through rings to find any five-membered ring that qualifies as a lactone
    for ring in ring_info:
        if len(ring) == 5:
            # Check each atom in the five-membered ring for the lactone features:
            # Identify a carbon that:
            #    1. Has a double bond to an oxygen that is not part of the ring (carbonyl oxygen)
            #    2. Is connected to at least one oxygen in the ring (ester oxygen)
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                # Look for a carbon atom candidate
                if atom.GetAtomicNum() == 6:
                    has_carbonyl = False
                    # Look through bonds of the candidate carbon
                    for bond in atom.GetBonds():
                        # Check for double bond to an oxygen atom
                        if bond.GetBondType() == Chem.BondType.DOUBLE:
                            neighbor = bond.GetOtherAtom(atom)
                            # The oxygen must be outside the ring (the carbonyl oxygen)
                            if neighbor.GetAtomicNum() == 8 and neighbor.GetIdx() not in ring:
                                has_carbonyl = True
                                break
                    
                    # If the candidate carbon does not have a proper carbonyl attachment, skip it.
                    if not has_carbonyl:
                        continue
                    
                    # Now verify that the candidate carbon is also directly bonded to an oxygen within the ring (ester oxygen)
                    has_ring_oxygen = False
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetIdx() in ring and neighbor.GetAtomicNum() == 8:
                            has_ring_oxygen = True
                            break
                    
                    if has_ring_oxygen and has_carbonyl:
                        return True, "Found five-membered lactone ring with ester connectivity (carbonyl and ring oxygen)."
    
    return False, "No gamma-lactone ring (five-membered lactone) found."

# Example usage (for testing purposes)
if __name__ == "__main__":
    # Example SMILES of a gamma-lactone: gamma-nonanolactone
    test_smiles = "CCCCCC[C@H]1CCC(=O)O1"
    result, reason = is_gamma_lactone(test_smiles)
    print("Result:", result)
    print("Reason:", reason)