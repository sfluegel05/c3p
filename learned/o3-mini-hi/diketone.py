"""
Classifies: CHEBI:46640 diketone
"""
"""
Classifies: A compound that contains two ketone functionalities (diketone).
A diketone here is defined as having exactly two ketone groups where each ketone group 
is a carbonyl group (C=O) whose carbon is bound (via single bonds) to two carbon atoms 
and has no attached hydrogen (so as to rule out aldehyde groups). Furthermore, we prefer to 
ignore carbonyl groups embedded in aromatic systems (e.g. quinones) as these often represent 
a different chemical class.
"""

from rdkit import Chem

def is_diketone(smiles: str):
    """
    Determines if a molecule is a diketone (exactly two independent ketone functionalities)
    based on its SMILES string. A ketone is defined here as a non-aromatic carbonyl carbon (C=O)
    that is connected via single bonds to two carbon atoms (i.e. it appears in an R-CO-R unit).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a diketone; False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # We'll search for ketone carbonyl carbons based on the following criteria:
    # 1. The atom must be a carbon ("C").
    # 2. It must have exactly one double-bond to oxygen.
    # 3. It must have exactly two other heavy-atom neighbors (which must be carbon)
    #    (so that the overall connectivity is: C(=O) with two C substituents).
    # 4. Its own atom environment should not be aromatic.
    #
    # Note: This heuristic may fail in special cases but should work for many acyclic and non‐quinone
    #       diketones.
    
    ketone_count = 0
    # To avoid double counting, we will count each qualifying carbonyl atom only once.
    for atom in mol.GetAtoms():
        # Only consider carbon atoms; skip if not "C"
        if atom.GetSymbol() != "C":
            continue
        # Exclude if the carbonyl carbon is aromatic.
        if atom.GetIsAromatic():
            continue
            
        # For a ketone carbon, there should be no attached hydrogen
        # (otherwise it would be an aldehyde, e.g. H-C=O).
        if atom.GetTotalNumHs() != 0:
            continue

        # Get all bonds around this carbon.
        bonds = atom.GetBonds()
        # We expect the carbonyl carbon (in a ketone: R-CO-R) to be connected to exactly three heavy atoms,
        # one from the C=O double bond and the other two being carbon atoms.
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        # If the number of heavy neighbors is not 3 then skip this atom.
        if len(heavy_neighbors) != 3:
            continue

        # Now check for exactly one double bond from this carbon to an oxygen.
        dbl_oxygen_found = False
        valid_ketone = True  # we assume valid unless one test fails
        # Count the number of bonds that are C=O double bonds.
        count_doubleO = 0
        for bond in bonds:
            # Only consider double bonds.
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                nbr = bond.GetOtherAtom(atom)
                if nbr.GetAtomicNum() == 8:
                    count_doubleO += 1
        if count_doubleO != 1:
            valid_ketone = False

        # Next, for all bonds that are not the C=O bond, they should connect only to carbon atoms.
        # (This ensures that we are not counting esters or acid derivatives.)
        if valid_ketone:
            for bond in bonds:
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    # Already accounted for C=O double bond.
                    continue
                nbr = bond.GetOtherAtom(atom)
                # We already ensured heavy_neighbors count above; now require that non-oxygen neighbors are carbon.
                if nbr.GetAtomicNum() != 6:
                    valid_ketone = False
                    break

        if valid_ketone:
            ketone_count += 1

    # Check whether the molecule has exactly two ketone functionalities.
    if ketone_count == 2:
        return True, "Contains exactly two ketone functionalities"
    else:
        return False, f"Found {ketone_count} ketone functionalities; exactly 2 are required for a diketone"

# Example usage; you can test with one of the provided SMILES.
if __name__ == "__main__":
    # Test with 2,3-Undecanedione: SMILES "O=C(CCCCCCCC)C(=O)C"
    test_smiles = "O=C(CCCCCCCC)C(=O)C"
    result, reason = is_diketone(test_smiles)
    print(f"Result: {result}\nReason: {reason}")