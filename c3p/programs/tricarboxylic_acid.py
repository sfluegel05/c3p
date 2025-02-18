"""
Classifies: CHEBI:27093 tricarboxylic acid
"""
"""
Classifies: Tricarboxylic Acid (An oxoacid containing three carboxy groups.)
"""
from rdkit import Chem

def is_tricarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a tricarboxylic acid based on its SMILES string.
    A tricarboxylic acid is defined as an oxoacid containing three carboxy groups.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a tricarboxylic acid, False otherwise
        str: Reason for the classification
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a carboxyl group.
    # This pattern covers both protonated and ionized carboxy groups.
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX1]")
    if carboxyl_pattern is None:
        return False, "Error creating SMARTS pattern for carboxyl group"
    
    # Find all substructure matches for the carboxyl group
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    n_carboxyl_groups = len(carboxyl_matches)
    
    # Check if exactly 3 carboxyl groups are found.
    if n_carboxyl_groups != 3:
        return False, f"Found {n_carboxyl_groups} carboxyl group(s); expected exactly 3"
    
    return True, "Contains exactly three carboxyl groups, consistent with a tricarboxylic acid"

# Example usage (uncomment below lines to test):
# smiles_examples = [
#     "OC(=O)CC(O)(CC(O)=O)C(O)=O",  # citric acid then some examples
#     "OC(=O)CCC(C(O)=O)C(=O)C(O)=O"  # another example
# ]
# for s in smiles_examples:
#     result, reason = is_tricarboxylic_acid(s)
#     print(s, result, reason)