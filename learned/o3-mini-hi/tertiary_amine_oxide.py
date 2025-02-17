"""
Classifies: CHEBI:134363 tertiary amine oxide
"""
#!/usr/bin/env python3
"""
Classifies: Tertiary amine oxide
Definition: An N-oxide in which a tetrahedral (sp3) nitrogen (with a +1 formal charge) is bonded to 
exactly one oxygen atom (with a -1 formal charge via a single bond) and three carbon atoms. Furthermore, 
each carbon substituent must be “simple” – that is, aside from its bond to the nitrogen, its immediate 
neighbors are only carbons or hydrogens (or, if in a ring, the entire ring must be free of heteroatoms).
"""

from rdkit import Chem

def is_tertiary_amine_oxide(smiles: str):
    """
    Determines if a molecule is a tertiary amine oxide based on its SMILES string.
    
    The improved algorithm:
      1. Parse the SMILES string.
      2. Iterate over atoms looking for nitrogen atoms that:
         - Have symbol "N" with a formal charge of +1.
         - Are tetrahedral (degree 4) and not aromatic.
      3. For each candidate nitrogen, require that exactly one neighbor is an oxygen with a -1 charge,
         that is single-bonded to N and that oxygen is linked only to that candidate nitrogen.
      4. The other three neighbors must be carbon atoms.
      5. For each carbon neighbor, check that its immediate neighbours (apart from the candidate nitrogen)
         consist only of carbons or hydrogens. If the carbon is in a ring, then if any other neighbor 
         “shares” a ring with the carbon the entire ring must include only carbons and hydrogens.
      6. If any candidate N passes these tests, return True with a success message.
         Otherwise, return False.

    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if a tertiary amine oxide center is found, False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information from the molecule (list of tuples of atom indices for each ring)
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # Helper function: given a candidate carbon (by its index) and one of its neighbor indices,
    # if that neighbor is not carbon or hydrogen then allow it only if the candidate carbon and that neighbor
    # appear together in a ring that contains exclusively carbon and hydrogen.
    def neighbor_allowed(candidate_idx, nbr):
        # direct allowed if the neighbor is carbon (atomic number 6) or hydrogen (atomic number 1)
        if nbr.GetAtomicNum() in (1, 6):
            return True

        # Otherwise, check if candidate and neighbor share a ring that is “simple”
        for ring in rings:
            if candidate_idx in ring and nbr.GetIdx() in ring:
                # Check if all atoms in this ring are carbons or hydrogens.
                if all(mol.GetAtomWithIdx(a).GetAtomicNum() in (1, 6) for a in ring):
                    return True
        return False

    # Iterate over candidate atoms
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "N" or atom.GetFormalCharge() != 1:
            continue
        # Require tetrahedral: degree equals 4 and not aromatic.
        if atom.GetDegree() != 4 or atom.GetIsAromatic():
            continue

        neighbors = atom.GetNeighbors()
        oxide_count = 0
        carbon_neighbors = []
        candidate_oxygen = None

        # Examine neighbors of the candidate nitrogen.
        for nbr in neighbors:
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            # Test if neighbor is an oxygen with proper formal charge and single-bonded.
            if nbr.GetSymbol() == "O" and nbr.GetFormalCharge() == -1:
                if bond is not None and bond.GetBondType() == Chem.BondType.SINGLE:
                    oxide_count += 1
                    candidate_oxygen = nbr
            elif nbr.GetAtomicNum() == 6:  # carbon atom
                carbon_neighbors.append(nbr)
            else:
                # Any other substituent type is not allowed under the strict definition.
                pass

        # We require exactly one oxygen and three carbon substituents.
        if oxide_count != 1 or len(carbon_neighbors) != 3:
            continue

        # Check the candidate oxygen: it should only be connected to this candidate nitrogen.
        if candidate_oxygen.GetDegree() != 1:
            continue

        # For each carbon neighbor, check that its immediate (1-bond) environment (ignoring the candidate N)
        # is only made of carbons or hydrogens (or if not, the offending neighbor is in a ring that is purely C/H).
        valid_substituents = True
        for c_atom in carbon_neighbors:
            for nbr in c_atom.GetNeighbors():
                if nbr.GetIdx() == atom.GetIdx():
                    continue
                if not neighbor_allowed(c_atom.GetIdx(), nbr):
                    valid_substituents = False
                    break
            if not valid_substituents:
                break

        if not valid_substituents:
            continue

        # If all tests pass for this candidate nitrogen then we classify the molecule as a tertiary amine oxide.
        return True, ("Found a tertiary amine oxide functional group "
                      "(N+ bonded to O- via a single bond and three carbon substituents, each of which is a "
                      "simple organic group)")

    return False, "No tertiary amine oxide functional group found"


# Optional testing routine when module is run directly.
if __name__ == "__main__":
    test_examples = [
        # True positives:
        ("C[N+](C)([O-])C", "trimethylamine N-oxide"),
        ("C1=CC([C@]2([N@+](CCC2)([O-])C)[H])=CN=C1", "trans-(S)-nicotine N(1')-oxide"),
        # False positive example:
        ("C(CC[N+](C)(C)[O-])N1C=2C(CCC3=C1C=CC=C3)=CC=CC2", "imipramine oxide (expected to be rejected)"),
    ]
    
    for s, name in test_examples:
        result, reason = is_tertiary_amine_oxide(s)
        print(f"SMILES: {s}\nName: {name}\nResult: {result}\nReason: {reason}\n")