"""
Classifies: CHEBI:82744 aliphatic aldoxime
"""
"""
Classifies: Aliphatic Aldoxime – any aldoxime derived from an aliphatic aldehyde.
A molecule is considered to have an aliphatic aldoxime group if it contains a functional group of 
the form R–CH=N–OH, where:
    - The carbon attached to the imine N has exactly one hydrogen (consistent with an aldehyde CH group).
    - That same carbon is not aromatic.
    - Its only heavy (non-hydrogen) neighbor (besides the imine nitrogen) is a carbon that is aliphatic (non‐aromatic).
"""

from rdkit import Chem

def is_aliphatic_aldoxime(smiles: str):
    """
    Determines if a molecule is an aliphatic aldoxime based on its SMILES string.
    Specifically, it checks for an oxime group of the form R–CH=N–OH where:
         - The carbon bonded to the imine nitrogen (CH) has exactly one hydrogen.
         - The candidate carbon is not aromatic.
         - Excluding the imine nitrogen and hydrogens, it has exactly one heavy neighbor,
           and that neighbor is an aliphatic (non-aromatic, atomic num 6) carbon.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as an aliphatic aldoxime, False otherwise.
        str: Detailed explanation of the decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to make hydrogen counting straightforward
    mol = Chem.AddHs(mol)
    
    # Define a SMARTS pattern to identify candidate oxime groups.
    # The pattern looks for a carbon double-bonded to a nitrogen that is connected to an -OH group: [CX3]=[NX2][OX2H]
    oxime_pattern = Chem.MolFromSmarts("[CX3]=[NX2][OX2H]")
    if oxime_pattern is None:
        return False, "Failed to compile SMARTS for oxime group"
    
    matches = mol.GetSubstructMatches(oxime_pattern)
    if not matches:
        return False, "No oxime group found"
    
    # Process each candidate match (order: carbon, nitrogen, oxygen)
    for match in matches:
        c_idx, n_idx, o_idx = match
        carbon = mol.GetAtomWithIdx(c_idx)
        nitrogen = mol.GetAtomWithIdx(n_idx)
        oxygen = mol.GetAtomWithIdx(o_idx)
        
        # Confirm that the oxygen is really part of a hydroxyl (-OH) (i.e. has at least one hydrogen)
        o_neighbors = oxygen.GetNeighbors()
        if not any(neigh.GetAtomicNum() == 1 for neigh in o_neighbors):
            continue  # Not a hydroxyl - skip

        # New check: ensure that the candidate carbon is not aromatic.
        if carbon.GetIsAromatic():
            continue
        
        # Count the hydrogen atoms attached to the candidate carbon.
        hydrogen_count = sum(1 for neigh in carbon.GetNeighbors() if neigh.GetAtomicNum() == 1)
        if hydrogen_count != 1:
            # For an aldehyde-derived oxime, the carbon should have exactly one hydrogen.
            continue
        
        # Examine the heavy (non-hydrogen) neighbors of the candidate carbon, excluding the imine nitrogen.
        heavy_neighbors = [neigh for neigh in carbon.GetNeighbors() 
                           if neigh.GetAtomicNum() != 1 and neigh.GetIdx() != n_idx]
        if len(heavy_neighbors) != 1:
            # There should be exactly one heavy substituent (the "R" group) other than the imine nitrogen.
            continue
        
        r_atom = heavy_neighbors[0]
        # The heavy neighbor should be a carbon atom.
        if r_atom.GetAtomicNum() != 6:
            continue
        # The heavy neighbor must be aliphatic (i.e. not aromatic)
        if r_atom.GetIsAromatic():
            continue
        
        # If the candidate passes all tests, then we have found an aliphatic aldoxime.
        return True, "Contains an aldoxime group (R–CH=N–OH) derived from an aliphatic aldehyde."
    
    # If no candidate oxime group qualifies, provide an explanation.
    return False, ("Oxime group(s) found but none appear to be derived from an aliphatic aldehyde "
                   "(the carbon does not have exactly one hydrogen, is aromatic, or its substituent is not an aliphatic carbon).")

# Example test cases (these can be commented out or modified as needed)
if __name__ == "__main__":
    test_smiles = [
        "C([C@@H](/C(=N\\O)/[H])C)C",          # (1Z,2S)-2-methylbutanal oxime
        "[H]\\C(C(C)C)=N/O",                   # (E)-2-methylpropanal oxime
        "C(\\CCCCCCCCSC)=N/O",                 # (E)-9-(methylsulfanyl)nonanal oxime
        "[H]C(=NO)C(C)CC",                     # 2-methylbutanal oxime
        "[H]\\C(=N/O)C(C)CC",                  # (1E,2S)-2-methylbutanal oxime
        "OC(C(O)C(O)\\C=N\\O)C(O)CO",           # (1E)-2,3,4,5,6-pentahydroxyhexanal oxime
        "[H]C(CCCCCCSC)=NO",                   # 7-(methylsulfanyl)heptanal oxime
        "C(CCCCCCCCSC)=NO",                    # 9-(methylsulfanyl)nonanal oxime
        "C(\\CCCCCCCSC)=N/O",                  # (E)-8-(methylsulfanyl)octanal oxime
        "[H]C(=NO)C(C)(C)SC",                  # 2-methyl-2-(methylsulfanyl)propanal oxime
        "C(CCCCSC)=NO",                        # 5-(methylsulfanyl)pentanal oxime
        "C(\\CCCCCSC)=N/O",                    # (E)-6-(methylsulfanyl)hexanal oxime
        "[H]C(CSC)=NO",                        # (methylsulfanyl)acetaldehyde oxime
        "[H]\\C(=N\\O)C(C)CC",                 # (Z)-2-methylbutanal oxime
        "[H]\\C(=N/O)C(C)CC",                  # (E)-2-methylbutanal oxime
        "C(CCCSC)=NO",                         # 4-(methylsulfanyl)butanal oxime
        "[H]\\C(C)=N\\O",                      # (Z)-acetaldehyde oxime
        "[H]\\C(C(C)C)=N\\O",                  # (Z)-2-methylpropanal oxime
        "[H]C(CCCCNC(=O)CCC(=O)N(O)CCCCCNC(=O)CCC(=O)N(O)CCCCCN)=NO",  # IC202C
        "C(\\CCCCSC)=N/O",                     # (E)-5-(methylsulfanyl)pentanal oxime
        "C(CCCCCSC)=NO",                       # 6-(methylsulfanyl)hexanal oxime
        "C(CCCCCCCSC)=NO",                     # 8-(methylsulfanyl)octanal oxime
        "[H]C(C)=NO"                           # acetaldehyde oxime
    ]
    
    for smi in test_smiles:
        result, reason = is_aliphatic_aldoxime(smi)
        print(f"SMILES: {smi}\nResult: {result}, Reason: {reason}\n")