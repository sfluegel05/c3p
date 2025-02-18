"""
Classifies: CHEBI:33704 alpha-amino acid
"""
"""
Classifies: CHEBI:33704 alpha-amino acid
"""
from rdkit import Chem

def is_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid based on its SMILES string.
    An alpha-amino acid has the amino group attached to the carbon adjacent to the carboxylic acid group.

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

    # Find carboxylic acid groups (-C(=O)OH or -C(=O)O-)
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches:
        return False, "No carboxylic acid group found"

    # Check each carboxylic acid group for adjacent amino group
    for carboxyl in carboxyl_matches:
        carboxylic_carbon = carboxyl[0]
        
        # Get neighboring atoms to carboxylic carbon
        neighbors = mol.GetAtomWithIdx(carboxylic_carbon).GetNeighbors()
        for neighbor in neighbors:
            # Check if neighbor is a carbon (alpha carbon)
            if neighbor.GetAtomicNum() == 6:
                alpha_carbon = neighbor
                # Check for amino group (-NH2, -NHR, -NR2) on alpha carbon
                alpha_neighbors = alpha_carbon.GetNeighbors()
                amino_groups = [n for n in alpha_neighbors if n.GetAtomicNum() == 7 and n.GetTotalNumHs() >= 1]
                
                if amino_groups:
                    return True, "Alpha carbon has amino group adjacent to carboxylic acid"

    return False, "No amino group found on alpha carbon relative to carboxylic acid"