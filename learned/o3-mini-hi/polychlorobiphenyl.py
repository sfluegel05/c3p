"""
Classifies: CHEBI:53156 polychlorobiphenyl
"""
"""
Classifies: Polychlorobiphenyl
Definition: A biphenyl compound containing between 2 and 10 chlorine atoms attached (exclusively) as substituents on the two benzene rings.
The improved approach identifies the biphenyl core and then verifies that every substituent on that core is chlorine.
"""

from rdkit import Chem

def is_polychlorobiphenyl(smiles: str):
    """
    Determines if a molecule is a polychlorobiphenyl based on its SMILES string.
    
    The function looks for a biphenyl scaffold built from two benzene rings connected by a single bond.
    For each candidate biphenyl core match, we examine the atoms directly attached to the core.
    Only chlorine substituents (atomic number 17) are allowed; any additional substituents lead to rejection.
    The candidate is accepted only if the number of chlorine substituents directly attached to the biphenyl core 
    is between 2 and 10 (inclusive).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule meets the criteria, False otherwise.
        str: Explanation for the classification decision.
    """
    
    # Construct the RDKit molecule from SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS for a biphenyl scaffold:
    # This pattern identifies two benzene rings (six-membered aromatic rings) connected by a single bond.
    biphenyl_smarts = "c1ccccc1-c2ccccc2"
    biphenyl_query = Chem.MolFromSmarts(biphenyl_smarts)
    
    # Search for biphenyl substructures in the molecule.
    matches = mol.GetSubstructMatches(biphenyl_query)
    if not matches:
        return False, "No biphenyl scaffold (two connected benzene rings) found"
    
    # Evaluate each biphenyl candidate match
    for match in matches:
        core_indices = set(match)  # set of atom indices involved in the biphenyl core (expected 12 atoms)
        chlorine_count = 0
        extra_substituent_found = False
        
        # Iterate over each atom in the biphenyl core.
        for idx in core_indices:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                # Only consider atoms not in the core (i.e. substituents)
                if nbr.GetIdx() not in core_indices:
                    # Check if this neighbor is chlorine (atomic number 17)
                    if nbr.GetAtomicNum() == 17:
                        chlorine_count += 1
                    else:
                        # Found a substituent that is not chlorine—reject this candidate.
                        extra_substituent_found = True
                        break
            if extra_substituent_found:
                break
        
        # If non-chlorine substituents are found, skip this candidate.
        if extra_substituent_found:
            continue
        
        # If only chlorine substituents were found, check if the chlorine count falls between 2 and 10.
        if 2 <= chlorine_count <= 10:
            return True, f"Contains a biphenyl scaffold with {chlorine_count} chlorine substituents on the core"
    
    # If none of the biphenyl match candidates pass the strict substituent test.
    return False, "No biphenyl core found with exclusively chlorine substituents on the rings (2 to 10 chlorines)"

# Example calls (for testing purposes):
if __name__ == "__main__":
    examples = [
        ("Clc1cc(Cl)c(Cl)c(c1)-c1cc(Cl)c(Cl)c1Cl", "2,2',3,3',5,5'-hexachlorobiphenyl (should be accepted)"),
        ("Clc1cc(-c2ccccc2)c(Cl)c(Cl)c1Cl", "2,3,4,5-tetrachlorobiphenyl (should be accepted)"),
        ("Oc1c(Cl)cc(Cl)cc1-c1ccccc1", "2-Hydroxy-3,5-dichlorobiphenyl (should be rejected due to extra substituent)"),
        ("Clc1ccc(cc1)-c1cc(Cl)c(Cl)c(Cl)c1", "2,3,4,3',4'-pentachlorobiphenyl (should be accepted)")
    ]
    for smi, desc in examples:
        result, reason = is_polychlorobiphenyl(smi)
        print(f"{desc}\n  SMILES: {smi}\n  Result: {result}\n  Reason: {reason}\n")