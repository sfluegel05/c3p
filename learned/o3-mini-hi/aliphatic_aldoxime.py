"""
Classifies: CHEBI:82744 aliphatic aldoxime
"""
"""
Classifies: Aliphatic Aldoxime – any aldoxime derived from an aliphatic aldehyde.
A molecule is considered to have an aliphatic aldoxime group if it contains a functional group of 
the form R–CH=N–OH, i.e. the carbon attached to =N must have exactly one hydrogen (CH) and one heavy
(aliphatic, non-aromatic carbon) substituent (R, besides the N).
The algorithm:
  1. Convert the SMILES to an RDKit molecule and add explicit hydrogens.
  2. Find candidate oxime groups by searching for a pattern where:
     - a carbon (sp2) forms a double bond to a nitrogen, and 
     - that nitrogen is single bonded to an oxygen which in turn bears a hydrogen (–OH).
  3. For each candidate the “aldehyde-derived” carbon is checked:
     - It must have exactly one hydrogen.
     - Among its heavy (non-hydrogen) neighbors, aside from the imine nitrogen, it must have exactly one substituent.
     - That substituent must be a carbon atom that is not aromatic.
  4. If any candidate passes these tests the function will return True; otherwise False.
"""

from rdkit import Chem

def is_aliphatic_aldoxime(smiles: str):
    """
    Determines if a molecule is an aliphatic aldoxime based on its SMILES string.
    The definition is: an aldoxime group of the form R–CH=N–OH, where:
         - the carbon (CH) has exactly one hydrogen (i.e. derived from an aldehyde)
         - and its R substituent (other neighbor of that carbon besides the imine-N) 
           is a heavy aliphatic (non-aromatic carbon).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as an aliphatic aldoxime, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse SMILES and add explicit hydrogens
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Define a SMARTS that finds candidate oxime groups.
    # This pattern matches a sp2 carbon (C) double-bonded to a nitrogen (N),
    # which is single-bonded to an oxygen (O) that in turn is attached to at least one hydrogen.
    # Our pattern: [CX3]=[NX2][OX2H]
    # Note: This pattern does not force a specific hydrogen count on the carbon.
    oxime_pattern = Chem.MolFromSmarts("[CX3]=[NX2][OX2H]")
    if oxime_pattern is None:
        return False, "Failed to compile SMARTS pattern for oxime group"
    
    # Find all matches to the candidate oxime pattern.
    # Each match is a tuple (c_idx, n_idx, o_idx) corresponding to the carbon, nitrogen, and oxygen.
    matches = mol.GetSubstructMatches(oxime_pattern)
    if not matches:
        return False, "No oxime group found"
    
    # Evaluate each candidate match.
    for match in matches:
        c_idx, n_idx, o_idx = match
        # Get the involved atoms.
        carbon = mol.GetAtomWithIdx(c_idx)
        nitrogen = mol.GetAtomWithIdx(n_idx)
        oxygen = mol.GetAtomWithIdx(o_idx)
        
        # Check if the oxygen is indeed -OH (i.e. bonded to at least one hydrogen).
        # (With AddHs the hydrogen(s) should be explicit.)
        o_neighbors = oxygen.GetNeighbors()
        has_hydrogen = any(neigh.GetAtomicNum() == 1 for neigh in o_neighbors)
        if not has_hydrogen:
            continue  # not a proper -OH group
        
        # For the aldehyde-derived carbon, count its hydrogens.
        # GetTotalNumHs returns the sum of implicit and explicit hydrogens (after AddHs all are explicit).
        num_h = carbon.GetTotalNumHs()
        if num_h != 1:
            # Some SMILES may include an extra explicitly written hydrogen (e.g. [H]C(...)=N/O)
            # in which case the optima structure after AddHs might show >1 hydrogen.
            # We require exactly one hydrogen for it to be derived from an aldehyde.
            continue
        
        # Now check the connectivity of the carbon.
        # It is double bonded to the nitrogen; aside from that, it should have one heavy substituent.
        heavy_neighbors = []
        for neigh in carbon.GetNeighbors():
            # Exclude the nitrogen that is part of the oxime group.
            if neigh.GetIdx() == n_idx:
                continue
            # Exclude hydrogens.
            if neigh.GetAtomicNum() == 1:
                continue
            heavy_neighbors.append(neigh)
        
        if len(heavy_neighbors) != 1:
            # If there are not exactly 1 heavy neighbor, then the underlying group is not a mono-substituted aldehyde.
            continue
        
        # Check that the substituent is an aliphatic carbon (atomic number 6) and not aromatic.
        r_group_atom = heavy_neighbors[0]
        if r_group_atom.GetAtomicNum() != 6:
            continue
        if r_group_atom.GetIsAromatic():
            continue
        
        # If this candidate passes all tests then it is an aliphatic aldoxime.
        return True, "Contains an aldoxime group derived from an aliphatic aldehyde (CH=N–OH with proper R-group)"
    
    # If we exhaust all candidate oxime groups without finding one that meets our criteria,
    # return a failure, providing a detailed reason.
    return False, ("Oxime group(s) found but none appear to be derived from an aliphatic aldehyde "
                   "(e.g., the carbon does not have exactly one hydrogen or its substituent is not an aliphatic carbon)")
    
    
# Example test cases (can be commented out if not needed)
if __name__ == "__main__":
    test_smiles = [
        "C([C@@H](/C(=N\\O)/[H])C)C",          # (1Z,2S)-2-methylbutanal oxime
        "[H]\\C(C(C)C)=N/O",                  # (E)-2-methylpropanal oxime
        "C(\\CCCCCCCCSC)=N/O",                # (E)-9-(methylsulfanyl)nonanal oxime
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
        "C(\\CCCCSC)=N/O",                    # (E)-5-(methylsulfanyl)pentanal oxime
        "C(CCCCCSC)=NO",                      # 6-(methylsulfanyl)hexanal oxime
        "C(CCCCCCCSC)=NO",                     # 8-(methylsulfanyl)octanal oxime
        "[H]C(C)=NO"                          # acetaldehyde oxime
    ]
    for smi in test_smiles:
        result, reason = is_aliphatic_aldoxime(smi)
        print(f"SMILES: {smi}\nResult: {result}, Reason: {reason}\n")