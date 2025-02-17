"""
Classifies: CHEBI:82744 aliphatic aldoxime
"""
"""
Classifies: Aliphatic Aldoxime – any aldoxime derived from an aliphatic aldehyde.
A molecule is considered an aliphatic aldoxime if it has an aldoxime group of the form R–CH=N–OH 
where the aldoxime carbon is CH (i.e. exactly one hydrogen) and the substituent R (the group attached
to the carbon aside from the oxime nitrogen) is a heavy (carbon) atom that is not aromatic.
"""

from rdkit import Chem

def is_aliphatic_aldoxime(smiles: str):
    """
    Determines if a molecule is an aliphatic aldoxime based on its SMILES string.
    The definition is: an aldoxime has a functional group of the form R–CH=N–OH (i.e. derived
    from an aliphatic aldehyde) where the aldehyde (CH) carbon bears exactly one hydrogen, and 
    the R-group (attached to that carbon aside from the N of the oxime) is a heavy aliphatic carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as an aliphatic aldoxime, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that we can count them directly.
    mol = Chem.AddHs(mol)
    
    # Define a SMARTS pattern for an aldoxime group derived from an aldehyde:
    # Using "[CH1](=N[OH])" forces the carbon to have exactly one hydrogen (CH).
    aldoxime_smarts = "[CH1](=N[OH])"
    aldoxime_pattern = Chem.MolFromSmarts(aldoxime_smarts)
    if aldoxime_pattern is None:
        return False, "Failed to compile SMARTS pattern for aldoxime"

    # Find all substructure matches for the aldoxime group.
    matches = mol.GetSubstructMatches(aldoxime_pattern)
    if not matches:
        return False, "No aldoxime group found"

    # Iterate over each aldoxime match. The match tuple contains indices: 
    # index 0: the aldehyde carbon, index 1: the imine nitrogen, index 2: the hydroxyl oxygen.
    for match in matches:
        aldehyde_carbon = mol.GetAtomWithIdx(match[0])
        # Although the SMARTS pattern "[CH1]" requires exactly one hydrogen by construction,
        # we double-check its total hydrogen count to be safe.
        if aldehyde_carbon.GetTotalNumHs() != 1:
            continue  # does not match an aldehyde-derived oxime

        # The carbon should not be aromatic.
        if aldehyde_carbon.GetIsAromatic():
            continue

        # From the aldehyde carbon, identify the neighboring heavy atom(s) excluding the oxime nitrogen.
        heavy_neighbors = []
        for neighbor in aldehyde_carbon.GetNeighbors():
            # Exclude the nitrogen that is part of the oxime (match[1])
            if neighbor.GetIdx() == match[1]:
                continue
            # Skip hydrogens
            if neighbor.GetAtomicNum() == 1:
                continue
            heavy_neighbors.append(neighbor)

        # For a typical aldehyde, the carbon should have exactly one heavy substituent (R group).
        if len(heavy_neighbors) != 1:
            continue

        r_neighbor = heavy_neighbors[0]
        # The substituent should be a carbon atom and should not be aromatic 
        # (i.e. the aldehyde comes from an aliphatic, not aromatic, aldehyde).
        if r_neighbor.GetAtomicNum() != 6:
            continue
        if r_neighbor.GetIsAromatic():
            continue

        # If we reach here, we have an aldoxime group with:
        # - a CH unit (exactly one hydrogen),
        # - the double bond to nitrogen and attached OH group,
        # - and the R substituent is a heavy aliphatic carbon.
        return True, "Contains an aldoxime group derived from an aliphatic aldehyde"

    # No valid aldoxime group met the criteria.
    return False, ("Aldoxime group(s) found but none appear to be derived from an aliphatic aldehyde "
                   "(e.g., aldehyde carbon does not have exactly one hydrogen or is attached to an aromatic substituent)")

# Example usage: (these test cases can be removed or commented out as needed)
if __name__ == "__main__":
    test_smiles = [
        "C([C@@H](/C(=N\\O)/[H])C)C",          # (1Z,2S)-2-methylbutanal oxime
        "[H]\\C(C(C)C)=N/O",                  # (E)-2-methylpropanal oxime
        "C(CCCCCCCCSC)=N/O",                  # (E)-9-(methylsulfanyl)nonanal oxime
        "[H]C(=NO)C(C)CC",                    # 2-methylbutanal oxime
        "[H]\\C(=N/O)C(C)CC",                 # (1E,2S)-2-methylbutanal oxime
        "OC(C(O)C(O)\\C=N\\O)C(O)CO",          # (1E)-2,3,4,5,6-pentahydroxyhexanal oxime
        "[H]C(CCCCCCSC)=NO",                  # 7-(methylsulfanyl)heptanal oxime
        "C(CCCCCCCCSC)=NO",                   # 9-(methylsulfanyl)nonanal oxime
        "C(\\CCCCCCCSC)=N/O",                 # (E)-8-(methylsulfanyl)octanal oxime
        "[H]C(=NO)C(C)(C)SC",                 # 2-methyl-2-(methylsulfanyl)propanal oxime
        "C(CCCCSC)=NO",                       # 5-(methylsulfanyl)pentanal oxime
        "C(\\CCCCCSC)=N/O",                   # (E)-6-(methylsulfanyl)hexanal oxime
        "[H]C(CSC)=NO",                       # (methylsulfanyl)acetaldehyde oxime
        "[H]\\C(=N\\O)C(C)CC",                # (Z)-2-methylbutanal oxime
        "[H]\\C(=N/O)C(C)CC",                 # (E)-2-methylbutanal oxime
        "C(CCCSC)=NO",                        # 4-(methylsulfanyl)butanal oxime
        "[H]\\C(C)=N\\O",                     # (Z)-acetaldehyde oxime
        "[H]\\C(C(C)C)=N\\O",                  # (Z)-2-methylpropanal oxime
        "[H]C(CCCCNC(=O)CCC(=O)N(O)CCCCCNC(=O)CCC(=O)N(O)CCCCCN)=NO",  # IC202C
    ]
    for smi in test_smiles:
        result, reason = is_aliphatic_aldoxime(smi)
        print(f"SMILES: {smi}\nResult: {result}, Reason: {reason}\n")