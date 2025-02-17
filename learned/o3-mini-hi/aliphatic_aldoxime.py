"""
Classifies: CHEBI:82744 aliphatic aldoxime
"""
"""
Classifies: Aliphatic Aldoxime – any aldoxime derived from an aliphatic aldehyde.
A molecule is considered to have an aliphatic aldoxime group if it contains a functional group of 
the form R–CH=N–OH where:
    - The carbon attached to the imine nitrogen has exactly one hydrogen (as in an aldehyde CH group).
    - That carbon is not aromatic.
    - Besides the imine nitrogen (and hydrogens), it has exactly one heavy neighbor.
    - The heavy substituent (the “R” group) is a carbon that is aliphatic (i.e. it and all connected heavy atoms are non‐aromatic).
"""

from rdkit import Chem

def is_aliphatic_aldoxime(smiles: str):
    """
    Determines if a molecule is an aliphatic aldoxime based on its SMILES string.
    It checks for an oxime group of the form R–CH=N–OH that comes from an aliphatic aldehyde.
    
    Requirements on the candidate aldehyde carbon (the one bonded to the imine N):
      - It is non-aromatic.
      - It has exactly one hydrogen.
      - Excluding the imine nitrogen and hydrogens, it has exactly one heavy (non-hydrogen) neighbor.
      - The attached heavy atom (which becomes the R group) must be a carbon and must be fully aliphatic;
        that is, by traversing all heavy atoms in its substituent, none of them should be aromatic.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an aliphatic aldoxime, False otherwise.
        str: Detailed explanation of the decision.
    """
    
    # Parse the SMILES string to create an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Adding explicit hydrogens to simplify hydrogen count checks.
    mol = Chem.AddHs(mol)
    
    # Define a SMARTS pattern to identify candidate oxime groups.
    # This pattern looks for a carbon double-bonded to a nitrogen which is bonded to an OH group.
    oxime_pattern = Chem.MolFromSmarts("[CX3]=[NX2][OX2H]")
    if oxime_pattern is None:
        return False, "Failed to compile SMARTS for oxime group"
    
    matches = mol.GetSubstructMatches(oxime_pattern)
    if not matches:
        return False, "No oxime group found"
    
    # Helper function: traverse the substituent (R-group) starting from the given heavy atom.
    # We exclude the candidate aldehyde carbon and all hydrogens.
    # If any heavy atom in this substituent is aromatic, we return False.
    def is_substituent_aliphatic(start_atom, exclude_idx):
        visited = set()
        stack = [start_atom]
        while stack:
            atom = stack.pop()
            if atom.GetIdx() in visited:
                continue
            visited.add(atom.GetIdx())
            if atom.GetIsAromatic():
                return False
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 1:
                    continue  # skip hydrogens
                if neighbor.GetIdx() == exclude_idx:
                    continue  # do not go back to the candidate aldehyde carbon
                if neighbor.GetIdx() not in visited:
                    stack.append(neighbor)
        return True

    # Process each candidate oxime match.
    # The pattern order is: candidate aldehyde carbon, imine nitrogen, and oxygen of OH.
    for match in matches:
        c_idx, n_idx, o_idx = match
        carbon = mol.GetAtomWithIdx(c_idx)
        nitrogen = mol.GetAtomWithIdx(n_idx)
        oxygen = mol.GetAtomWithIdx(o_idx)
        
        # Verify the oxygen is truly part of a hydroxyl (-OH) group (i.e., it has at least one hydrogen).
        if not any(neigh.GetAtomicNum() == 1 for neigh in oxygen.GetNeighbors()):
            continue  # skip candidate
        
        # The candidate aldehyde carbon must be non-aromatic.
        if carbon.GetIsAromatic():
            continue
        
        # Count attached hydrogens for the candidate carbon.
        hydrogen_count = sum(1 for neigh in carbon.GetNeighbors() if neigh.GetAtomicNum() == 1)
        if hydrogen_count != 1:
            # For an aldehyde-derived oxime, the carbon should have exactly one hydrogen.
            continue
        
        # Examine heavy (non-hydrogen) neighbors excluding the imine nitrogen.
        heavy_neighbors = [neigh for neigh in carbon.GetNeighbors() 
                           if neigh.GetAtomicNum() != 1 and neigh.GetIdx() != n_idx]
        if len(heavy_neighbors) != 1:
            # There should be exactly one heavy substituent (the R group).
            continue
        
        r_atom = heavy_neighbors[0]
        # The heavy neighbor should be a carbon.
        if r_atom.GetAtomicNum() != 6:
            continue
        
        # Check that the immediate substituent carbon is non-aromatic.
        if r_atom.GetIsAromatic():
            continue
        
        # Traverse the substituent: if any heavy atom is aromatic, reject this candidate.
        if not is_substituent_aliphatic(r_atom, exclude_idx=carbon.GetIdx()):
            continue
        
        # If all conditions are met then we have found an aliphatic aldoxime.
        return True, "Contains an aldoxime group (R–CH=N–OH) derived from an aliphatic aldehyde."
        
    # No candidate oxime group qualifies after the additional aliphatic check.
    return False, ("Oxime group(s) found but none appear to be derived from an aliphatic aldehyde "
                   "(the candidate carbon or its substituent does not meet aliphatic criteria).")

# Example test cases
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