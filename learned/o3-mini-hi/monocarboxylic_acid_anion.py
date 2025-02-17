"""
Classifies: CHEBI:35757 monocarboxylic acid anion
"""
#!/usr/bin/env python3
"""
Classifies: Monocarboxylic acid anion
Definition: A carboxylic acid anion formed when the carboxy group of a monocarboxylic acid is deprotonated.
The molecule should contain exactly one deprotonated carboxyl group, represented as "[CX3](=O)[O-]".
"""

from rdkit import Chem

def is_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a monocarboxylic acid anion based on its SMILES string.
    
    A monocarboxylic acid anion is defined as a molecule forming from a monocarboxylic acid upon the
    loss of its acidic proton – in structure, it should contain exactly one carboxylate group, i.e., a 
    carbonyl (C=O) directly bonded to an anionic oxygen [O-].
    
    Args:
        smiles (str): The SMILES string representing the molecule.
        
    Returns:
        (bool, str): A tuple where the first element is True if the molecule is a monocarboxylic acid anion,
                     False otherwise, and the second element is a string that provides the reasoning for the result.
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for a deprotonated carboxylic acid group.
    # [CX3] ensures that the carboxyl carbon is sp2 hybridized.
    # The pattern matches a carbonyl group (C=O) directly bonded to an oxygen with a negative charge.
    carboxylate_smarts = "[CX3](=O)[O-]"
    carboxylate_pattern = Chem.MolFromSmarts(carboxylate_smarts)
    if carboxylate_pattern is None:
        return False, "Error creating SMARTS pattern for carboxylate group"
        
    # Find all matches in the molecule
    matches = mol.GetSubstructMatches(carboxylate_pattern)
    
    if len(matches) == 0:
        return False, "No deprotonated carboxyl (carboxylate) group found"
    elif len(matches) > 1:
        return False, f"Found {len(matches)} carboxylate groups; molecule is not a monocarboxylic acid anion"
    
    # Additional sanity check: ensure that the deprotonated group is not part of a larger, polyacidic scaffold.
    # We can check the total formal charge on the molecule. Typically, for a monocarboxylate anion, the overall charge should be -1.
    total_charge = sum([atom.GetFormalCharge() for atom in mol.GetAtoms()])
    if total_charge != -1:
        return False, f"Expected overall charge of -1 for a monocarboxylate anion, found charge = {total_charge}"
    
    return True, "The molecule contains exactly one deprotonated carboxyl group and has an overall charge of -1, consistent with a monocarboxylic acid anion"

# Example usage (you can remove these lines if using as module)
if __name__ == "__main__":
    test_smiles = "[O-]C(=O)c1nc2[nH]c(=O)[nH]c(=O)c2[nH]1"  # xanthine-8-carboxylate example
    result, reason = is_monocarboxylic_acid_anion(test_smiles)
    print(result, reason)