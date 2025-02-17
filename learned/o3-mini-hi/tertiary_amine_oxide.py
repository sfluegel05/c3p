"""
Classifies: CHEBI:134363 tertiary amine oxide
"""
#!/usr/bin/env python3
"""
Classifies: Tertiary amine oxide
Definition: An N-oxide in which a tetrahedral (sp³) nitrogen (with a +1 formal charge) 
is bonded to exactly one oxygen atom (with a -1 formal charge via a single bond) and three carbon atoms.
Furthermore, each carbon substituent must be “simple” – aside from its bond to the nitrogen, 
its immediate neighbors must consist solely of carbons or hydrogens. If a substituent neighbor is not carbon or hydrogen, 
it is acceptable only if the bond comes as part of a ring that is exclusively composed of carbons and hydrogens.
"""

from rdkit import Chem

def is_tertiary_amine_oxide(smiles: str):
    """
    Determines if a molecule is a tertiary amine oxide based on its SMILES string.
    
    Improved Algorithm:
      1. Parse the SMILES string.
      2. Iterate over all atoms looking for nitrogen atoms that satisfy:
         - Symbol "N" with a formal charge of +1.
         - Tetrahedral (degree of 4) and not aromatic.
      3. For each candidate N, count the neighbors:
         - Exactly one neighbor must be oxygen having a -1 charge, connected by a SINGLE bond.
         - The oxygen must be connected only to this candidate (degree=1).
         - The other three neighbors must be carbon atoms.
      4. For each carbon neighbor, check its immediate neighborhood (excluding the candidate nitrogen):
         - Every neighbor must be either carbon (atomic number 6) or hydrogen (atomic number 1).
         - If any neighbor is not C or H, then check if the substituent carbon and that neighbor share a ring that is composed exclusively of C and H.
         - If no such ring exists for any “non-simple” neighbor, the substituent is disqualified.
      5. If at least one candidate nitrogen passes all tests, return True with an explanation.
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
    
    # Retrieve ring information as a list of tuples containing atom indices per ring.
    rings = mol.GetRingInfo().AtomRings()

    # Improved helper function: Check if a neighbor (nbr) of a candidate carbon (by its index) is allowed.
    # Allowed if the neighbor is C or H or if there exists at least one ring that contains both the candidate carbon
    # and this neighbor and that ring is exclusively composed of C and H.
    def neighbor_allowed(carbon_idx, nbr):
        # If the neighbor is carbon or hydrogen, it is simple.
        if nbr.GetAtomicNum() in (1, 6):
            return True
        # Otherwise, check for a ring shared between the candidate carbon and the neighbor.
        for ring in rings:
            if carbon_idx in ring and nbr.GetIdx() in ring:
                # Ring is allowed only if every atom in it is either C or H.
                if all(mol.GetAtomWithIdx(a).GetAtomicNum() in (1,6) for a in ring):
                    return True
        return False

    # Iterate over all atoms looking for a candidate nitrogen.
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "N" or atom.GetFormalCharge() != 1:
            continue
        # Require tetrahedral geometry: a degree of 4 and not aromatic.
        if atom.GetDegree() != 4 or atom.GetIsAromatic():
            continue

        neighbors = atom.GetNeighbors()
        oxide_count = 0
        carbon_neighbors = []
        candidate_oxygen = None

        # Analyze neighbors of the candidate nitrogen.
        for nbr in neighbors:
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            # Check if neighbor is oxygen with formal charge -1 and is single-bonded.
            if nbr.GetSymbol() == "O" and nbr.GetFormalCharge() == -1:
                if bond is not None and bond.GetBondType() == Chem.BondType.SINGLE:
                    oxide_count += 1
                    candidate_oxygen = nbr
            elif nbr.GetAtomicNum() == 6:  # Must be a carbon atom.
                carbon_neighbors.append(nbr)
            else:
                # Any substituent that is not oxygen or carbon is not acceptable.
                pass

        # The candidate N must have exactly one oxygen and three carbon substituents.
        if oxide_count != 1 or len(carbon_neighbors) != 3:
            continue

        # Check that the candidate oxygen is bonded only to this nitrogen.
        if candidate_oxygen.GetDegree() != 1:
            continue

        # For each carbon neighbor, ensure its immediate surroundings (excluding the N) are simple.
        valid_substituents = True
        for c_atom in carbon_neighbors:
            for nbr in c_atom.GetNeighbors():
                if nbr.GetIdx() == atom.GetIdx():
                    continue  # Skip back-reference to candidate nitrogen.
                # If a neighbor is not allowed, mark the substituent as complex.
                if not neighbor_allowed(c_atom.GetIdx(), nbr):
                    valid_substituents = False
                    break
            if not valid_substituents:
                break

        if not valid_substituents:
            continue  # This candidate nitrogen did not pass the substituent simplicity check.

        # If we reached here, we found a tertiary amine oxide center.
        return True, ("Found a tertiary amine oxide functional group "
                      "(N+ bonded to O- via a single bond and three carbon substituents, each of which is a "
                      "simple organic group)")

    # No candidate nitrogen satisfied all criteria.
    return False, "No tertiary amine oxide functional group found"


# Optional: Only run tests if this module is executed directly.
if __name__ == "__main__":
    test_examples = [
        # True positives (from provided examples):
        ("[H][C@]1(C[N+]2([O-])CC[C@@]34C(Nc5ccccc35)=C(C(=O)OC)[C@@]1([H])C[C@]24[H])C(C)=O", "(-)-Alstolucine C"),
        ("C[N+](C)([O-])C", "trimethylamine N-oxide"),
        ("C1=CC([C@]2([N@+](CCC2)([O-])C)[H])=CN=C1", "trans-(S)-nicotine N(1')-oxide"),
        # False positive example:
        ("C(CC[N+](C)(C)[O-])N1C=2C(CCC3=C1C=CC=C3)=CC=CC2", "imipramine oxide (should not be classified)"),
    ]
    
    for s, name in test_examples:
        result, explanation = is_tertiary_amine_oxide(s)
        print(f"SMILES: {s}\nName: {name}\nResult: {result}\nExplanation: {explanation}\n")