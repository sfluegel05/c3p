"""
Classifies: CHEBI:134363 tertiary amine oxide
"""
#!/usr/bin/env python3
"""
Classifies: Tertiary amine oxide
Definition: An N-oxide where there are three organic groups bonded to the nitrogen atom.
A tertiary amine oxide is characterized by a nitrogen atom with a positive formal charge,
bonded to an oxygen atom with a negative formal charge, and three additional (organic) groups.
"""

from rdkit import Chem

def is_tertiary_amine_oxide(smiles: str):
    """
    Determines if a molecule is a tertiary amine oxide based on its SMILES string.

    Steps:
      1. Parse the molecule from SMILES.
      2. Iterate over all nitrogen atoms.
      3. For each nitrogen with a +1 formal charge that is bonded to exactly four atoms,
         check if exactly one neighbor is an oxygen with a -1 formal charge.
      4. Verify that the other three substituents are organic groups (i.e., not hydrogens).
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule contains at least one tertiary amine oxide center, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate over atoms in the molecule to find a candidate nitrogen center
    for atom in mol.GetAtoms():
        # Look for nitrogen with a +1 formal charge.
        if atom.GetSymbol() == "N" and atom.GetFormalCharge() == 1:
            neighbors = list(atom.GetNeighbors())
            # For a tertiary amine oxide, nitrogen should be tetravalent
            if len(neighbors) != 4:
                continue

            oxide_count = 0
            organic_group_count = 0
            for nbr in neighbors:
                # Check if neighbor is the oxide group: oxygen with -1 charge.
                if nbr.GetSymbol() == "O" and nbr.GetFormalCharge() == -1:
                    oxide_count += 1
                else:
                    # Ensure that the substituent is not just a hydrogen.
                    if nbr.GetAtomicNum() != 1:
                        organic_group_count += 1
            # The functional group should have exactly one oxide and three organic substituents.
            if oxide_count == 1 and organic_group_count == 3:
                return True, "Found a tertiary amine oxide functional group (N+ bonded to O- and three organic groups)"
    
    return False, "No tertiary amine oxide functional group found"

# Example usage (for testing purposes):
if __name__ == "__main__":
    test_smiles = [
        "C[N+](C)([O-])C",  # trimethylamine N-oxide
        "CN(C)(=O)c1ccccc1",  # NOT a tertiary amine oxide (this is an N-oxide in an amide environment)
    ]
    for s in test_smiles:
        result, reason = is_tertiary_amine_oxide(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")