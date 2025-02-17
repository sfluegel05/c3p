"""
Classifies: CHEBI:82744 aliphatic aldoxime
"""
"""
Classifies: Aliphatic Aldoxime – any aldoxime derived from an aliphatic aldehyde.
"""

from rdkit import Chem

def is_aliphatic_aldoxime(smiles: str):
    """
    Determines if a molecule is an aliphatic aldoxime based on its SMILES string.
    An aliphatic aldoxime has an aldoxime group (RCH=N–OH) where the aldehyde carbon 
    (R–CH) bears exactly one hydrogen and is not part of an aromatic system.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic aldoxime, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to ensure reliable hydrogen counting.
    mol = Chem.AddHs(mol)

    # Define a SMARTS pattern for an aldoxime functional group:
    # [CX3](=N[OX2H]) corresponds to a carbon (from the aldehyde) double-bonded to N,
    # which in turn is bonded to an -OH group.
    aldoxime_smarts = "[CX3](=N[OX2H])"
    aldoxime_pattern = Chem.MolFromSmarts(aldoxime_smarts)
    if aldoxime_pattern is None:
        return False, "Failed to compile SMARTS pattern for aldoxime"

    # Find all substructure matches of the aldoxime pattern
    matches = mol.GetSubstructMatches(aldoxime_pattern)
    if not matches:
        return False, "No aldoxime group found"

    # Iterate over found aldoxime groups. The SMARTS pattern has 3 atoms:
    # index 0: the carbon (aldehyde carbon), index 1: the imine nitrogen, index 2: the hydroxyl oxygen.
    valid_match_found = False
    for match in matches:
        aldehyde_carbon = mol.GetAtomWithIdx(match[0])
        # The aldehyde carbon should have exactly one hydrogen (as in RCH=NOH)
        if aldehyde_carbon.GetTotalNumHs() != 1:
            continue  # not an aldehyde-derived oxime

        # The aldehyde carbon must be aliphatic (not in an aromatic system)
        if aldehyde_carbon.GetIsAromatic():
            continue

        # Identify the substituent R (the neighbor that is not the nitrogen in the oxime group)
        r_neighbor = None
        for neighbor in aldehyde_carbon.GetNeighbors():
            # Skip hydrogen atoms and the nitrogen atom that is part of the oxime
            if neighbor.GetAtomicNum() == 1:
                continue
            if neighbor.GetIdx() == match[1]:
                continue
            r_neighbor = neighbor
            break
        
        if r_neighbor is None:
            # If we do not find a heavy-atom neighbor aside from the oxime nitrogen,
            # then the structure does not follow R–CH=N–OH.
            continue

        # Check the nature of the substituent: if it is a carbon atom, it should not be aromatic.
        if r_neighbor.GetAtomicNum() == 6 and r_neighbor.GetIsAromatic():
            continue

        # If we passed all tests, we have found an aldoxime group derived from an aliphatic aldehyde.
        valid_match_found = True
        break

    if not valid_match_found:
        return False, ("Aldoxime group(s) found but none appear to be derived from an aliphatic aldehyde "
                       "(e.g., aldehyde carbon does not have exactly one hydrogen or is attached to an aromatic ring)")

    return True, "Contains an aldoxime group derived from an aliphatic aldehyde"

# Example usage (you can remove or comment these out when integrating into larger projects):
if __name__ == "__main__":
    # Testing with some of the provided examples.
    test_smiles = [
        "C([C@@H](/C(=N\\O)/[H])C)C",          # (1Z,2S)-2-methylbutanal oxime
        "[H]C(C)=NO",                         # acetaldehyde oxime
        "[H]C(CCCCNC(=O)CCC(=O)N(O)CCCCCNC(=O)CCC(=O)N(O)CCCCCN)=NO"  # IC202C
    ]
    for smi in test_smiles:
        result, reason = is_aliphatic_aldoxime(smi)
        print(f"SMILES: {smi}\nResult: {result}, Reason: {reason}\n")