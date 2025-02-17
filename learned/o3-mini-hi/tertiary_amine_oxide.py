"""
Classifies: CHEBI:134363 tertiary amine oxide
"""
#!/usr/bin/env python3
"""
Classifies: Tertiary amine oxide
Definition: An N-oxide where there are three organic groups (i.e., carbon‐based substituents)
bonded to the nitrogen atom. In our interpretation the nitrogen should be saturated (sp3)
with a formal +1 charge; it is bonded to exactly one oxygen atom (having a –1 charge and linked
by a single bond) and three carbon atoms.
"""

from rdkit import Chem

def is_tertiary_amine_oxide(smiles: str):
    """
    Determines if a molecule is a tertiary amine oxide based on its SMILES string.
    
    The algorithm:
      1. Parse the SMILES into an RDKit molecule.
      2. For each nitrogen atom in the molecule:
          a. Require it has a formal charge of +1 and exactly four bonds (making it tetrahedral).
          b. Exactly one neighbor must be an oxygen with formal charge –1 and the N–O bond must be single.
          c. The other three neighbors must be carbon atoms (organic substituents).
      3. If any such nitrogen atom is found, return True along with a success message.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if a tertiary amine oxide center is found, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse SMILES into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Iterate over all atoms looking for candidate nitrogen centers
    for atom in mol.GetAtoms():
        # Look for nitrogen atoms with +1 charge.
        if atom.GetSymbol() != "N" or atom.GetFormalCharge() != 1:
            continue
        # For a tertiary amine oxide, the N should be tetrahedral (4 neighbors) 
        # and not aromatic (typically sp3 rather than sp2)
        if atom.GetDegree() != 4 or atom.GetIsAromatic():
            continue
        
        neighbors = atom.GetNeighbors()
        oxide_count = 0
        organic_count = 0
        
        # Check each neighbor
        for nbr in neighbors:
            # Get the bond between the candidate N and this neighbor.
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            if nbr.GetSymbol() == "O" and nbr.GetFormalCharge() == -1:
                # Require that the bond is single
                if bond is not None and bond.GetBondType() == Chem.BondType.SINGLE:
                    oxide_count += 1
            # For organic groups, we require that the neighbor is a carbon atom.
            elif nbr.GetAtomicNum() == 6:
                organic_count += 1
            # Otherwise, if the substituent is not a carbon and not the oxide,
            # then it does not count as an organic substituent.
        
        # Check if we have exactly one oxide neighbor and three carbons
        if oxide_count == 1 and organic_count == 3:
            return True, "Found a tertiary amine oxide functional group (N+ bonded to O- via a single bond and three carbon substituents)"
    
    return False, "No tertiary amine oxide functional group found"

# Example usage for quick testing (this will only run when the module is executed directly)
if __name__ == "__main__":
    test_examples = [
        # True positives:
        "C[N+](C)([O-])C",  # trimethylamine N-oxide
        "C1=CC([C@]2([N@+](CCC2)([O-])C)[H])=CN=C1",  # trans-(S)-nicotine N(1')-oxide
        # False positive example (should be rejected by our stricter requirements):
        "C(CC[N+](C)(C)[O-])N1C=2C(CCC3=C1C=CC=C3)=CC=CC2"  # imipramine oxide (expected to be rejected)
    ]
    
    for s in test_examples:
        result, reason = is_tertiary_amine_oxide(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")