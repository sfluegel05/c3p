"""
Classifies: CHEBI:53156 polychlorobiphenyl
"""
"""
Classifies: Polychlorobiphenyl
Definition: A biphenyl compound containing between 2 and 10 chlorine atoms attached to the two benzene rings.
The new approach identifies the biphenyl core using a SMARTS search and then counts the chlorine atoms directly bonded to the core.
Any extra substituents (non-chlorine) found on the core atoms are tolerated.
"""

from rdkit import Chem

def is_polychlorobiphenyl(smiles: str):
    """
    Determines if a molecule is a polychlorobiphenyl based on its SMILES string.
    
    The function looks for a biphenyl scaffold built from two benzene rings connected by a single bond.
    Then it counts the number of chlorine atoms directly attached to the atoms in the biphenyl core.
    If the chlorine count is between 2 and 10 (inclusive), the molecule is classified as a polychlorobiphenyl.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule meets the criteria, False otherwise.
        str: Explanation for the classification decision.
    """
    
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS for a biphenyl scaffold:
    # It matches two benzene rings (six-membered aromatic rings) connected by a single bond.
    biphenyl_smarts = "c1ccccc1-c2ccccc2"
    biphenyl_query = Chem.MolFromSmarts(biphenyl_smarts)
    
    # Search for the biphenyl substructure in the molecule.
    matches = mol.GetSubstructMatches(biphenyl_query)
    if not matches:
        return False, "No biphenyl scaffold (two connected benzene rings) found"
    
    # For each match candidate, get the core atom indices (first 6 belong to ring1 and next 6 to ring2)
    for match in matches:
        core_indices = set(match)  # set of indices for the biphenyl core atoms
        chlorine_count = 0
        
        # Count chlorine atoms (atomic number 17) directly bonded to any core atom.
        # We iterate over each core atom and then check its neighbors
        for idx in core_indices:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                # Only consider neighbors not in the core.
                if nbr.GetIdx() not in core_indices:
                    if nbr.GetAtomicNum() == 17:  # chlorine
                        chlorine_count += 1
        # Check if the chlorine count is between 2 and 10.
        if 2 <= chlorine_count <= 10:
            return True, f"Contains a biphenyl scaffold with {chlorine_count} chlorine substituents on the core"
    
    # If none of the biphenyl core matches has the appropriate chlorine substitution, return False.
    return False, "No biphenyl core found with the required chlorine substitution pattern (2 to 10 chlorine atoms)"

# Example calls (for testing purposes):
if __name__ == "__main__":
    examples = [
        ("Clc1cc(Cl)c(Cl)c(c1)-c1cc(Cl)c(Cl)c1Cl", "2,2',3,3',5,5'-hexachlorobiphenyl (should be accepted)"),
        ("Clc1cc(-c2ccccc2)c(Cl)c(Cl)c1Cl", "2,3,4,5-tetrachlorobiphenyl (should be accepted)"),
        ("Clc1cccc(c1)-c1cc(Cl)c(Cl)c(Cl)c1Cl", "2,3,3',4,5-pentachlorobiphenyl (should be accepted)"),
        ("Clc1ccc(cc1Cl)-c1ccc(Cl)c(Cl)c1Cl", "2,3,4,3',4'-pentachlorobiphenyl (should be accepted)"),
        ("OC1=C(C=C(C=C1[N+]([O-])=O)Cl)C2=C(O)C([N+](=O)[O-])=CC(=C2)Cl", "niclofolan (should be accepted if the core meets the criteria)")
    ]
    for smi, desc in examples:
        result, reason = is_polychlorobiphenyl(smi)
        print(f"{desc}\n  SMILES: {smi}\n  Result: {result}\n  Reason: {reason}\n")