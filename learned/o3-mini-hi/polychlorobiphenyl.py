"""
Classifies: CHEBI:53156 polychlorobiphenyl
"""
"""
Classifies: Polychlorobiphenyl
Definition: A biphenyl compound containing between 2 and 10 chlorine atoms attached to the two benzene rings.
This improved approach first identifies a biphenyl scaffold (two benzene rings connected by a single bond),
and then counts the chlorine atoms attached directly to the ring atoms (i.e. not in the core).
Any additional substituents do not lead to rejection. The accepted candidate must have the chlorine count (on the rings)
between 2 and 10 (inclusive).
"""

from rdkit import Chem

def is_polychlorobiphenyl(smiles: str):
    """
    Determines if a molecule is a polychlorobiphenyl based on its SMILES string.
    
    The function looks for a biphenyl scaffold in the molecule. It then counts the number of chlorine atoms (atomic number 17)
    that are directly bonded to the atoms of the biphenyl core. The molecule is accepted if at least one candidate biphenyl core
    has between 2 and 10 chlorine substituents.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets the criteria, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string to create an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS for a biphenyl scaffold:
    # This pattern finds two benzene rings (aromatic hexagons) connected by a single bond.
    biphenyl_smarts = "c1ccccc1-c2ccccc2"
    biphenyl_query = Chem.MolFromSmarts(biphenyl_smarts)
    
    # Search for biphenyl matches in the molecule.
    matches = mol.GetSubstructMatches(biphenyl_query)
    if not matches:
        return False, "No biphenyl scaffold (two connected benzene rings) found"
    
    # For each biphenyl match candidate, count the chlorine atoms attached to the core.
    for match in matches:
        core_indices = set(match)  # indices of the atoms in the biphenyl core (expected 12 atoms)
        chlorine_count = 0
        
        # For each atom in the biphenyl core, look at neighbors that are not in the core.
        for idx in core_indices:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() not in core_indices:
                    # Count this substituent if it is chlorine
                    if nbr.GetAtomicNum() == 17:
                        chlorine_count += 1
        # Check if the number of chlorine substituents is between 2 and 10.
        if 2 <= chlorine_count <= 10:
            return True, f"Contains a biphenyl scaffold with {chlorine_count} chlorine substituents on the rings"
    
    # If none of the biphenyl matches had the appropriate number of chlorine substituents.
    return False, "No biphenyl core found with 2 to 10 chlorine substituents attached to the rings"

# Example calls for testing purposes:
if __name__ == "__main__":
    examples = [
        ("Clc1cc(Cl)c(Cl)c(c1)-c1cc(Cl)c(Cl)c1Cl", "2,2',3,3',5,5'-hexachlorobiphenyl (expected: True)"),
        ("Clc1cc(-c2ccccc2)c(Cl)c(Cl)c1Cl", "2,3,4,5-tetrachlorobiphenyl (expected: True)"),
        ("Clc1cccc(c1)-c1cc(Cl)c(Cl)c(Cl)c1Cl", "2,3,3',4,5-pentachlorobiphenyl (expected: True)"),
        ("ClC1=C(OC2=C(Cl)C=C(Cl)C=C2)C(C3=C(O)C(Cl)=CC(=C3)Cl)=CC(=C1O)Cl", 
         "Ambigol E (expected: True even with extra substituents)"),
        ("OC1=C(C=C(C=C1[N+]([O-])=O)Cl)C2=C(O)C([N+](=O)[O-])=CC(=C2)Cl",
         "niclofolan (expected: True if biphenyl core is identified with 2-10 chlorines)")
    ]
    for smi, desc in examples:
        result, reason = is_polychlorobiphenyl(smi)
        print(f"{desc}\n  SMILES: {smi}\n  Result: {result}\n  Reason: {reason}\n")