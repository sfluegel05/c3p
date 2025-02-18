"""
Classifies: CHEBI:53156 polychlorobiphenyl
"""
"""
Classifies: Polychlorobiphenyl
Definition: A biphenyl compound containing between 2 and 10 chlorine atoms attached to the two benzene rings.
Additional Requirement: The only substituents on the benzene rings must be chlorine atoms (aside from hydrogens).
"""

from rdkit import Chem

def is_polychlorobiphenyl(smiles: str):
    """
    Determines if a molecule is a polychlorobiphenyl based on its SMILES string.
    
    Requirements for a polychlorobiphenyl:
      - It must contain a biphenyl scaffold: two connected benzene rings.
      - On the benzene rings of the scaffold, the only substituents (apart from hydrogens) allowed are chlorine atoms.
      - The total number of chlorine atoms attached to the biphenyl core is between 2 and 10.
    
    Args:
       smiles (str): SMILES string of the molecule.
       
    Returns:
       bool: True if the molecule meets the criteria, False otherwise.
       str: Explanation for the classification decision.
    """
    # Parse SMILES and add explicit hydrogens so our substituents are clearly visible.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Define a SMARTS pattern for a biphenyl scaffold (two benzene rings connected by a single bond).
    biphenyl_pattern = Chem.MolFromSmarts("c1ccc(cc1)-c2ccccc2")
    matches = mol.GetSubstructMatches(biphenyl_pattern)
    if not matches:
        return False, "Biphenyl scaffold not detected"
    
    # Use the first substructure match to define the biphenyl core.
    # (The match returns a tuple of atom indices belonging to the two benzene rings.)
    biphenyl_atoms = set(matches[0])
    
    # Now inspect each aromatic carbon in the biphenyl core:
    #   a) Count chlorine atoms attached directly to the benzene carbons.
    #   b) Ensure that no substituent (other than hydrogen or chlorine) is attached.
    chlorine_count = 0
    for idx in biphenyl_atoms:
        atom = mol.GetAtomWithIdx(idx)
        # Check that the atom is aromatic and a carbon.
        if not (atom.GetAtomicNum() == 6 and atom.GetIsAromatic()):
            continue
        # Now examine each neighbor that is not part of the biphenyl core.
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() in biphenyl_atoms:
                continue
            # Allowed substituents: hydrogen or chlorine.
            atomic_num = neighbor.GetAtomicNum()
            if atomic_num == 17:
                chlorine_count += 1
            elif atomic_num != 1:
                # Found a substituent other than Cl or H on the biphenyl core.
                return False, f"Disallowed substituent ({neighbor.GetSymbol()}) found on the biphenyl core"
    
    if chlorine_count < 2:
        return False, f"Found {chlorine_count} chlorine atoms attached to the biphenyl scaffold; at least 2 are required"
    if chlorine_count > 10:
        return False, f"Found {chlorine_count} chlorine atoms attached to the biphenyl scaffold; no more than 10 are allowed"
    
    return True, f"Contains a biphenyl scaffold with {chlorine_count} chlorine substituents and no extra groups on the rings"

# Example calls (for testing purpose):
if __name__ == "__main__":
    examples = [
        ("Clc1cc(Cl)c(Cl)c(c1)-c1cc(Cl)cc(Cl)c1Cl", "2,2',3,3',5,5'-hexachlorobiphenyl"),
        ("Oc1c(Cl)cc(Cl)cc1-c1ccccc1", "2-Hydroxy-3,5-dichlorobiphenyl (should be rejected)"),
    ]
    for smi, name in examples:
        result, reason = is_polychlorobiphenyl(smi)
        print(f"{name} -> {result}: {reason}")