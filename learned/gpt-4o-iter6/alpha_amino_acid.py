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
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches:
        return False, "No carboxyl group found"

    # Look for an alpha-carbon bonded to both carboxyl (C) and amino group (N)
    for carboxyl_match in carboxyl_matches:
        carboxyl_c = carboxyl_match[0]  # C in 'C(=O)O'
        # Get neighbors of this C to find the alpha carbon
        neighbors = [nbr.GetIdx() for nbr in mol.GetAtomWithIdx(carboxyl_c).GetNeighbors() if nbr.GetAtomicNum() == 6]
        
        for alpha_c in neighbors:
            # Get neighbors of the alpha carbon and check if any is a nitrogen
            alpha_neighbors = [nbr.GetIdx() for nbr in mol.GetAtomWithIdx(alpha_c).GetNeighbors()]
            if any(mol.GetAtomWithIdx(nbr).GetAtomicNum() == 7 for nbr in alpha_neighbors):  # N is atomic number 7
                return True, "Contains carboxyl group and amino group at alpha position"

    return False, "No amino group at alpha position found"