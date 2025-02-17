"""
Classifies: CHEBI:32877 primary amine
"""
"""
Classifies: Primary amine
A primary amine is defined as a compound formally derived from ammonia (NH3)
by replacing one hydrogen atom with a hydrocarbyl group (i.e. R–NH2).
"""

from rdkit import Chem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule is a primary amine based on its SMILES string.
    A primary amine should contain at least one nitrogen atom that is bonded
    to exactly one non-hydrogen (preferably a carbon) substituent and two hydrogens.
    Amide nitrogen atoms (bonded to C=O) are not considered primary amines.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a primary amine group, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Iterate over all atoms in the molecule to search for primary amine groups.
    for atom in mol.GetAtoms():
        # Check if the atom is nitrogen (atomic number 7)
        if atom.GetAtomicNum() != 7:
            continue
        
        # Get the total number of hydrogens (implicit + explicit)
        num_H = atom.GetTotalNumHs()
        
        # Get heavy atom neighbors (non-hydrogen atoms)
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        
        # For a primary amine group (R-NH2):
        # The nitrogen should have exactly one heavy neighbor and two hydrogens.
        if len(heavy_neighbors) != 1 or num_H != 2:
            continue
        
        # Check that the substituent attached is a carbon (the hydrocarbyl group)
        substituent = heavy_neighbors[0]
        if substituent.GetAtomicNum() != 6:
            continue
        
        # Exclude cases where the nitrogen is part of an amide (e.g. adjacent to a C=O group)
        is_amide = False
        # Look at bonds from the substituent carbon
        for bond in substituent.GetBonds():
            # Only consider double bonds
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                # Get the atom at the other end of the bond from the substituent carbon
                other = bond.GetOtherAtom(substituent)
                if other.GetAtomicNum() == 8:
                    # Found a C=O, likely an amide if nitrogen is attached
                    is_amide = True
                    break
        if is_amide:
            continue

        # If we reach here, the nitrogen qualifies as a primary amine (R-NH2)
        return True, "Contains a primary amine group (R-NH2)"
    
    # If no qualifying nitrogen is found, return False.
    return False, "No primary amine group (R-NH2) found"

# Example usage and testing on a few provided SMILES strings:
if __name__ == "__main__":
    test_smiles = [
        "Nc1ccccc1",                        # aniline, primary aromatic amine
        "C[C@@H](N)Cc1ccccc1",               # (R)-amphetamine
        "NCCCCCCCCCCCCCC",                  # tetradecylamine (long-chain primary amine)
        "CN"                                # methylamine, primary amine
    ]
    
    for smi in test_smiles:
        result, reason = is_primary_amine(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")