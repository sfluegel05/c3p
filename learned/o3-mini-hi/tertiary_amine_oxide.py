"""
Classifies: CHEBI:134363 tertiary amine oxide
"""
#!/usr/bin/env python3
"""
Classifies: Tertiary amine oxide
Definition: An N-oxide where there are three organic groups (i.e., carbon‐based substituents)
bonded to the nitrogen atom. In our interpretation the nitrogen should be tetrahedral (sp3) with
a formal +1 charge; it must be bonded to exactly one oxygen atom (having a –1 charge and connected
by a single bond) and three carbon atoms. In addition, each carbon substituent should be “simple”
in that its immediate neighbors (apart from the amine-oxide nitrogen) are only carbons or hydrogens.
This extra check helps eliminate candidates embedded in complex heterocyclic or heteroatom–rich 
environments.
"""

from rdkit import Chem

def is_tertiary_amine_oxide(smiles: str):
    """
    Determines if a molecule is a tertiary amine oxide based on its SMILES string.
    
    The improved algorithm:
      1. Parse the SMILES into an RDKit molecule.
      2. Iterate over all atoms looking for nitrogen atoms that:
          a. Have symbol "N" with a formal charge of +1.
          b. Are tetrahedral (degree 4) and not aromatic.
      3. For each such candidate nitrogen, count the neighbors.
         - Exactly one neighbor must be an oxygen atom with a formal charge of –1,
           and the N–O bond must be single; also require that the oxygen is only bonded
           to this nitrogen.
         - The other three neighbors must be carbon atoms.
      4. Further, each carbon neighbor is examined: in its own 1-bond environment (excluding 
         the candidate nitrogen), all attached heavy atoms must be carbons (atomic number 6). 
         This is a simple heuristic to enforce that the substituent is a “pure” organic group.
      5. If any candidate nitrogen satisfies these conditions, return True with a success message.
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

    # Iterate over all atoms looking for candidate nitrogen centers
    for atom in mol.GetAtoms():
        # Check if the atom is nitrogen with formal charge +1
        if atom.GetSymbol() != "N" or atom.GetFormalCharge() != 1:
            continue
        # Nitrogen should have exactly 4 neighbors (tetrahedral) and not be aromatic.
        if atom.GetDegree() != 4 or atom.GetIsAromatic():
            continue

        neighbors = atom.GetNeighbors()
        oxide_count = 0
        organic_count = 0

        # To store the candidate oxygen and carbon neighbors
        candidate_oxygen = None
        candidate_carbons = []
        
        for nbr in neighbors:
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            # Check for oxygen neighbor with formal charge -1 and single bond
            if nbr.GetSymbol() == "O" and nbr.GetFormalCharge() == -1:
                if bond is not None and bond.GetBondType() == Chem.BondType.SINGLE:
                    oxide_count += 1
                    candidate_oxygen = nbr
            elif nbr.GetAtomicNum() == 6:  # Carbon atom
                organic_count += 1
                candidate_carbons.append(nbr)
            # Any other neighbor is not allowed for a tertiary amine oxide center.

        # Basic check for counts: exactly one oxide and three carbons.
        if oxide_count != 1 or organic_count != 3:
            continue

        # Additional check: ensure the oxygen is connected only to the candidate nitrogen.
        if candidate_oxygen is not None:
            # Get neighbors of the oxygen (should only be the N).
            o_neighbors = candidate_oxygen.GetNeighbors()
            if len(o_neighbors) != 1:
                continue

        # Further check: for each carbon neighbor, verify that its other neighbors (except the candidate N)
        # consist only of carbons or hydrogens.
        valid_substituents = True
        for c_atom in candidate_carbons:
            for nbr in c_atom.GetNeighbors():
                # Skip the candidate nitrogen itself
                if nbr.GetIdx() == atom.GetIdx():
                    continue
                # Allow hydrogen (atomic number 1) and carbon (atomic number 6)
                if nbr.GetAtomicNum() not in (1, 6):
                    valid_substituents = False
                    break
            if not valid_substituents:
                break

        if not valid_substituents:
            continue

        # If all tests pass for this candidate N, we classify the molecule as a tertiary amine oxide.
        return True, ("Found a tertiary amine oxide functional group "
                      "(N+ bonded to O- via a single bond and three carbon substituents "
                      "with each carbon substituent being a simple organic group)")
    
    return False, "No tertiary amine oxide functional group found"

# Optional: quick testing when running the module directly.
if __name__ == "__main__":
    test_examples = [
        # True positives:
        "C[N+](C)([O-])C",  # trimethylamine N-oxide
        "C1=CC([C@]2([N@+](CCC2)([O-])C)[H])=CN=C1",  # trans-(S)-nicotine N(1')-oxide
        # False positive example:
        "C(CC[N+](C)(C)[O-])N1C=2C(CCC3=C1C=CC=C3)=CC=CC2"  # imipramine oxide (expected to be rejected)
    ]
    
    for s in test_examples:
        result, reason = is_tertiary_amine_oxide(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")