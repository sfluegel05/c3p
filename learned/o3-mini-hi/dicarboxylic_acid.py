"""
Classifies: CHEBI:35692 dicarboxylic acid
"""
"""
Classifies: Any carboxylic acid containing two carboxy groups (dicarboxylic acid)
"""

from rdkit import Chem

def is_dicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a dicarboxylic acid based on its SMILES string.
    A dicarboxylic acid is defined as any carboxylic acid containing exactly two carboxy groups.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a dicarboxylic acid, False otherwise.
        str: Explanation or reason for the classification decision.
    """
    
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a carboxylic acid group.
    # This pattern corresponds to the –C(=O)–OH group.
    carboxyl_smarts = "[CX3](=O)[OX2H]"
    carboxyl_pattern = Chem.MolFromSmarts(carboxyl_smarts)
    if carboxyl_pattern is None:
        return False, "Error in generating SMARTS pattern for carboxylic acid"

    # Find all non-overlapping substructure matches for the carboxyl group
    matches = mol.GetSubstructMatches(carboxyl_pattern)
    num_groups = len(matches)
    
    # Classification: must contain exactly two carboxy groups.
    if num_groups == 2:
        return True, "Molecule contains exactly two carboxylic acid groups"
    elif num_groups < 2:
        return False, f"Found only {num_groups} carboxylic acid group(s); need exactly two."
    else:
        return False, f"Found {num_groups} carboxylic acid groups; dicarboxylic acid must have exactly two."
        
# Example testing
if __name__ == "__main__":
    # Example: meso-tartaric acid (SMILES: O[C@@H]([C@@H](O)C(O)=O)C(O)=O)
    test_smiles = "O[C@@H]([C@@H](O)C(O)=O)C(O)=O"
    result, reason = is_dicarboxylic_acid(test_smiles)
    print(f"SMILES: {test_smiles}\nResult: {result}\nReason: {reason}")