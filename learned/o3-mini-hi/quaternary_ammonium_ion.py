"""
Classifies: CHEBI:35267 quaternary ammonium ion
"""
"""
Classifies: Quaternary ammonium ion
Definition: A derivative of ammonium (NH4+) in which all four hydrogens bonded to nitrogen 
have been replaced with univalent (usually organyl) groups.
Our approach:
- Parse the SMILES string to an RDKit molecule.
- Look for any nitrogen atom that has a formal charge of +1, a total bonding count (degree) of exactly 4,
  and no hydrogen atoms attached (explicit or implicit).
- Additionally, for our improved screening we require that all four substituents attached to the nitrogen are carbon-based.
  This helps to ensure that the positive nitrogen comes from substitution of ammonium (NH4+) rather than some other motif.
  
Note: This approach may be too strict for unusual quaternary ammonium ions, but it aims to reduce the false positives observed.
"""

from rdkit import Chem

def is_quaternary_ammonium_ion(smiles: str):
    """
    Checks if a molecule contains at least one quaternary ammonium ion.
    
    In this improved code, the charged nitrogen must:
      - have a formal charge of +1,
      - be bonded to exactly 4 atoms (no implicit hydrogens),
      - have no hydrogen atoms attached (i.e. fully substituted),
      - each substituent must be carbon (atomic number 6), corresponding to organyl groups.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if a quaternary ammonium ion is found, False otherwise.
        str: A reason explaining the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Iterate over all atoms in the molecule
    for atom in mol.GetAtoms():
        # We are looking for a nitrogen atom (atomic num 7) with a positive formal charge (+1)
        if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1:
            # In a proper quaternary ammonium ion, the nitrogen should have 4 connections, and no hydrogen atoms.
            # RDKit assigns implicit hydrogens if necessary so GetTotalNumHs() should be zero here.
            if atom.GetDegree() == 4 and atom.GetTotalNumHs() == 0:
                # Now check that every neighbor is a carbon; this is our attempt to enforce “organyl” substituents.
                neighbors = atom.GetNeighbors()
                # If any neighbor is not a carbon (atomic number 6), then skip this atom.
                if all(neighbor.GetAtomicNum() == 6 for neighbor in neighbors):
                    return True, "Found a nitrogen atom with +1 charge, four carbon-based substituents, and no hydrogens (quaternary ammonium ion)"
    
    return False, "No quaternary ammonium group (N+ with 4 carbon-based substituents) found"

# Example usage:
if __name__ == "__main__":
    # A test case known to contain a quaternary ammonium group:
    test_smiles = "C[N+](C)(C)CC(O)O"  # betaine aldehyde hydrate
    result, reason = is_quaternary_ammonium_ion(test_smiles)
    print("Result:", result)
    print("Reason:", reason)