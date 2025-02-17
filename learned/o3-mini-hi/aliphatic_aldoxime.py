"""
Classifies: CHEBI:82744 aliphatic aldoxime
"""
"""
Classifies: Aliphatic Aldoxime – any aldoxime derived from an aliphatic aldehyde.
A molecule is considered to have an aliphatic aldoxime group if it contains a functional group of 
the form R–CH=N–OH, where the carbon (CH) has exactly one hydrogen (consistent with an aldehyde),
and its other neighbor (besides the imine N) is an aliphatic (non-aromatic carbon).
"""

from rdkit import Chem

def is_aliphatic_aldoxime(smiles: str):
    """
    Determines if a molecule is an aliphatic aldoxime based on its SMILES string.
    The specific functionality searched is an aldoxime group of the form R–CH=N–OH, in which:
         - The carbon attached to the imine N must have exactly one hydrogen.
         - Excluding the imine N and any hydrogens, the only other substituent should be
           an aliphatic carbon (i.e. atomic number 6 and not aromatic).
    
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
    
    # Define a SMARTS pattern to find candidate oxime groups.
    # We look for a carbon (sp2) double-bonded to a nitrogen, which in turn is
    # single-bonded to an oxygen that has at least one hydrogen: [CX3]=[NX2][OX2H]
    oxime_pattern = Chem.MolFromSmarts("[CX3]=[NX2][OX2H]")
    if oxime_pattern is None:
        return False, "Failed to compile SMARTS for oxime group"
    
    matches = mol.GetSubstructMatches(oxime_pattern)
    if not matches:
        return False, "No oxime group found"
    
    # Process each candidate match
    for match in matches:
        # Match ordering: (carbon, nitrogen, oxygen)
        c_idx, n_idx, o_idx = match
        carbon = mol.GetAtomWithIdx(c_idx)
        nitrogen = mol.GetAtomWithIdx(n_idx)
        oxygen = mol.GetAtomWithIdx(o_idx)
        
        # Confirm that oxygen is part of a hydroxyl group (-OH)
        o_neighbors = oxygen.GetNeighbors()
        if not any(neigh.GetAtomicNum() == 1 for neigh in o_neighbors):
            continue  # If oxygen isn't attached to any hydrogen, skip
        
        # Instead of relying on GetTotalNumHs (which can be skewed by explicit H in SMILES),
        # we count the explicit hydrogen atoms bonded to our candidate carbon.
        hydrogen_count = sum(1 for neigh in carbon.GetNeighbors() if neigh.GetAtomicNum() == 1)
        if hydrogen_count != 1:
            # There must be exactly one hydrogen attached to the candidate carbon
            continue
        
        # Now examine the neighbors of the candidate carbon.
        # It is already bonded to the imine nitrogen and one hydrogen.
        # If we exclude the nitrogen (n_idx) and hydrogens, we should have exactly one heavy neighbor,
        # and that neighbor should be a carbon (atomic num 6) and not aromatic.
        heavy_neighbors = [neigh for neigh in carbon.GetNeighbors() 
                           if neigh.GetAtomicNum() != 1 and neigh.GetIdx() != n_idx]
        if len(heavy_neighbors) != 1:
            continue
        
        r_atom = heavy_neighbors[0]
        if r_atom.GetAtomicNum() != 6:
            continue
        if r_atom.GetIsAromatic():
            continue
            
        # If the candidate passes all tests, we have found an aliphatic aldoxime.
        return True, "Contains an aldoxime group (R–CH=N–OH) derived from an aliphatic aldehyde."
    
    # If no candidate oxime group qualifies, provide an explanation.
    return False, ("Oxime group(s) found but none appear to be derived from an aliphatic aldehyde "
                   "(the carbon does not have exactly one hydrogen or its substituent is not a non-aromatic carbon).")
                   
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
        "[H]\\C(C(C)C)=N\\O",                   # (Z)-2-methylpropanal oxime
        "[H]C(CCCCNC(=O)CCC(=O)N(O)CCCCCNC(=O)CCC(=O)N(O)CCCCCN)=NO",  # IC202C
        "C(\\CCCCSC)=N/O",                     # (E)-5-(methylsulfanyl)pentanal oxime
        "C(CCCCCSC)=NO",                       # 6-(methylsulfanyl)hexanal oxime
        "C(CCCCCCCSC)=NO",                      # 8-(methylsulfanyl)octanal oxime
        "[H]C(C)=NO"                           # acetaldehyde oxime
    ]
    
    for smi in test_smiles:
        result, reason = is_aliphatic_aldoxime(smi)
        print(f"SMILES: {smi}\nResult: {result}, Reason: {reason}\n")