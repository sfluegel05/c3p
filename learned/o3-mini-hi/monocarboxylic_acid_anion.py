"""
Classifies: CHEBI:35757 monocarboxylic acid anion
"""
#!/usr/bin/env python3
"""
Classifies: Monocarboxylic acid anion
Definition: A carboxylic acid anion formed when the carboxyl group of a monocarboxylic acid is deprotonated.
That is, the molecule should have exactly one deprotonated carboxyl group, leading to the overall -1 charge,
and no other charged atoms outside of that carboxylate group.
"""

from rdkit import Chem

def is_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a monocarboxylic acid anion based on its SMILES string.
    
    A monocarboxylic acid anion should have exactly one carboxylate group (defined as [CX3](=O)[O-])
    and an overall formal charge of -1. Additionally, to avoid misclassifications, none of the other atoms
    in the molecule should have a nonzero formal charge.
    
    Args:
        smiles (str): The SMILES string representing the molecule.
        
    Returns:
        (bool, str): A tuple where the first element is True if the molecule is classified as a 
                     monocarboxylic acid anion, and False otherwise. The second element provides
                     an explanation for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for a deprotonated carboxyl group: the carboxylate.
    carboxylate_smarts = "[CX3](=O)[O-]"
    carboxylate_pattern = Chem.MolFromSmarts(carboxylate_smarts)
    if carboxylate_pattern is None:
        return False, "Error creating SMARTS pattern for carboxylate group"
    
    # Find all matches of the carboxylate pattern.
    matches = mol.GetSubstructMatches(carboxylate_pattern)
    if len(matches) == 0:
        return False, "No deprotonated carboxyl (carboxylate) group found"
    elif len(matches) > 1:
        return False, f"Found {len(matches)} carboxylate groups; molecule is not a monocarboxylic acid anion"
    
    # Check overall formal charge
    total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if total_charge != -1:
        return False, f"Expected overall charge of -1 for a monocarboxylate anion, found charge = {total_charge}"
    
    # Get the indices of atoms in the matched carboxylate group.
    carboxylate_atoms = set(matches[0])
    
    # Now check that no other atom (outside the carboxylate group) carries any formal charge.
    # This helps eliminate zwitterions or molecules with other ionized groups.
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() != 0 and atom.GetIdx() not in carboxylate_atoms:
            return False, ("Additional charged atoms found outside the deprotonated carboxyl group, "
                           "indicating the presence of other ionizable groups")
    
    return True, ("The molecule contains exactly one deprotonated carboxyl group and has an overall charge of -1 "
                  "with no other charged atoms, consistent with a monocarboxylic acid anion")

# Example usage (can be removed if used as a module)
if __name__ == "__main__":
    test_smiles = "[O-]C(=O)c1nc2[nH]c(=O)[nH]c(=O)c2[nH]1"  # xanthine-8-carboxylate example
    result, reason = is_monocarboxylic_acid_anion(test_smiles)
    print(result, reason)