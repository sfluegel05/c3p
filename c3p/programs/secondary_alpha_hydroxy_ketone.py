"""
Classifies: CHEBI:2468 secondary alpha-hydroxy ketone
"""
"""
Classifies: Secondary Alpha-Hydroxy Ketone (Acyloin)

Definition:
    A secondary α‐hydroxy ketone (acyloin) is a compound containing the motif 
        R–CH(OH)–C(=O)–R′ 
    in which the CH bearing –OH is secondary (bonded to three heavy atoms: one from –OH, one from C=O, and one organyl group).
    
This implementation parses the molecule and then iterates over carbon atoms to verify:
  1. The candidate carbon has exactly three heavy atom neighbors.
  2. Exactly one neighbor is an -OH group (an oxygen carrying at least one hydrogen),
  3. One neighbor is a carbonyl carbon (a carbon double bonded to an oxygen),
  4. The remaining neighbor is a carbon (the organyl group).
  5. The candidate carbon has exactly one hydrogen (explicit+implicit).
  
If a candidate meeting these criteria is found, the molecule is classified as containing the secondary alpha-hydroxy ketone motif.
"""

from rdkit import Chem

def is_secondary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is a secondary alpha-hydroxy ketone (acyloin) based on its SMILES string.
    The motif must be: R–CH(OH)–C(=O)–R′, where the CH bearing -OH is secondary (three heavy neighbors).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule contains the secondary alpha-hydroxy ketone motif, False otherwise.
        str : Reason for classification.
    """
    # Parse SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Iterate over candidate atoms; we look for a carbon that might be the acyloin center.
    for atom in mol.GetAtoms():
        # Check if the atom is carbon.
        if atom.GetAtomicNum() != 6:
            continue
        
        # Check that the carbon has exactly one hydrogen (explicit+implicit)
        if atom.GetTotalNumHs() != 1:
            continue
        
        # Count heavy atom neighbors (neighbors with atomic number > 1).
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        if len(heavy_neighbors) != 3:
            continue
        
        # Initialize counters for the three neighbor roles:
        oh_found = False
        carbonyl_found = False
        organyl_found = False

        # We iterate over the neighbors.
        for nbr in atom.GetNeighbors():
            # Case 1: Check for hydroxyl group: oxygen attached that in turn has at least one hydrogen.
            if nbr.GetAtomicNum() == 8:
                # We require this oxygen not to be part of a carbonyl (i.e. not double-bonded to a carbon)
                # and to have at least one hydrogen.
                # Get total hydrogen count on this oxygen.
                if nbr.GetTotalNumHs() >= 1:
                    oh_found = True
                    continue

            # Case 2: Check for carbonyl group: a carbon neighbor that has a double bond to oxygen.
            if nbr.GetAtomicNum() == 6:
                # Check if this carbon (nbr) has a double bond to an oxygen.
                found_dbl_oxygen = False
                for bond in nbr.GetBonds():
                    # Get the neighbor on the other side of the bond.
                    other = bond.GetOtherAtom(nbr)
                    if other.GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() == 2.0:
                        found_dbl_oxygen = True
                        break
                if found_dbl_oxygen:
                    carbonyl_found = True
                    continue
            
            # Case 3: The remaining neighbor should be an organic substituent, typically a carbon.
            if nbr.GetAtomicNum() == 6:
                organyl_found = True
                continue

        # After checking all three neighbors, confirm that we found all required groups.
        if oh_found and carbonyl_found and organyl_found:
            return True, "Molecule contains the secondary alpha-hydroxy ketone (acyloin) motif"
    
    return False, "Does not contain the required secondary alpha-hydroxy ketone motif"

# Optional: execute tests when run as a script.
if __name__ == '__main__':
    # Test example: (S)-benzoin should be classified as an acyloin.
    test_smiles = "O[C@H](C(=O)c1ccccc1)c1ccccc1"
    result, reason = is_secondary_alpha_hydroxy_ketone(test_smiles)
    print(result, reason)