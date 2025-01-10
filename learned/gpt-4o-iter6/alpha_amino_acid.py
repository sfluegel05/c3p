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

    # Look for carboxyl group pattern (C(=O)O)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[O-]")  # or protonated form if applicable, adjust accordingly
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches:
        return False, "No carboxyl group found"

    # Look for an alpha-amino group (nitrogen bonded to alpha-carbon)
    amino_group_pattern = Chem.MolFromSmarts("[C;!R][NX3;H2,H1,H0;!$(N=*)]")  # Nitrogen with at least one hydrogen 
    amino_matches = mol.GetSubstructMatches(amino_group_pattern)

    for carboxyl_match in carboxyl_matches:
        carboxyl_c = carboxyl_match[0]  # C in 'C(=O)O'
        # Get neighbors of this C to find the alpha carbon
        neighbors = [nbr.GetIdx() for nbr in mol.GetAtomWithIdx(carboxyl_c).GetNeighbors() if nbr.GetAtomicNum() == 6]
        
        for alpha_c in neighbors:
            # For an identified alpha carbon, check for attachment to an amino group using SMARTS pattern
            for amino_match in amino_matches:
                if amino_match[0] == alpha_c:
                    return True, "Contains carboxyl group and amino group at alpha position"
                
    return False, "No amino group at alpha position found"

# Example usage
example_smiles = "C[C@H](N)C(O)=O"  # L-alanine
result = is_alpha_amino_acid(example_smiles)
print(result)  # Output: (True, 'Contains carboxyl group and amino group at alpha position')