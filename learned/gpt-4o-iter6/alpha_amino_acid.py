"""
Classifies: CHEBI:33704 alpha-amino acid
"""
from rdkit import Chem

def is_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid based on its SMILES string.
    An alpha-amino acid is defined as an amino acid in which the amino group 
    is located on the carbon atom at the position alpha to the carboxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxy group pattern (C(=O)O)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxyl group found"

    # Look for alpha-amino group pattern (N attached to the carbon next to carboxyl group)
    # Define the patterns for fragments:
    # C[N] - an amino group must be attached directly to a carbon
    c_n_pattern = Chem.MolFromSmarts("[C][NX3][H]") # Nitrogen with three connections and at least one hydrogen (NH2)
    
    # The amino group must be alpha to the carboxy group
    for carboxyl_match in mol.GetSubstructMatches(carboxyl_pattern):      
        carboxyl_c = carboxyl_match[0] # C in 'C(=O)O'

        # Get neighbors of this C, find the alpha carbon
        neighbors = mol.GetAtomWithIdx(carboxyl_c).GetNeighbors()
        for neighbor in neighbors:
            if neighbor.GetAtomicNum() == 6:  # Check if it is carbon (alpha-carbon)
                # Check it has a nitrogen attached
                alpha_c = neighbor.GetIdx()
                alpha_c_neighbors = neighbor.GetNeighbors()
                for nn in alpha_c_neighbors:
                    if nn.GetAtomicNum() == 7 and mol.HasSubstructMatch(c_n_pattern, nn.GetIdx()): 
                        return True, "Contains carboxyl group and amino group at alpha position"
                
    return False, "No amino group at alpha position found"

# Example usage
example_smiles = "C[C@H](N)C(O)=O"  # L-alanine
result = is_alpha_amino_acid(example_smiles)
print(result)  # Output: (True, 'Contains carboxyl group and amino group at alpha position')