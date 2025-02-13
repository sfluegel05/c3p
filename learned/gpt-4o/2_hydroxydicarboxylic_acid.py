"""
Classifies: CHEBI:50263 2-hydroxydicarboxylic acid
"""
from rdkit import Chem

def is_2_hydroxydicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxydicarboxylic acid based on its SMILES string.
    A 2-hydroxydicarboxylic acid is a dicarboxylic acid with a hydroxy group on the alpha carbon to one of the carboxy groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 2-hydroxydicarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string to an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the carboxylic acid and hydroxy groups
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    hydroxy_pattern = Chem.MolFromSmarts("O")

    # Identify all carboxylic acids in the molecule
    carboxylic_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_matches) < 2:
        return False, f"Found {len(carboxylic_matches)} carboxylic acid groups, need at least 2"

    # Check for hydroxy group on the alpha carbon to a carboxylic group
    for carbonylic_group in carboxylic_matches:
        # Get the carbon atom index in carboxylic acid
        carboxylic_carbon_idx = carbonylic_group[0]
        
        # Check neighbors of this carbon for hydroxy pattern
        for neighbor in mol.GetAtomWithIdx(carboxylic_carbon_idx).GetNeighbors():
            if neighbor.GetSymbol() == 'C':
                # Check if this neighbor carbon has a hydroxy group attached
                alpha_carbon_idx = neighbor.GetIdx()
                if mol.GetAtomWithIdx(alpha_carbon_idx).GetTotalNumHs():
                    if mol.GetSubstructMatches(Chem.MolFromSmarts(f"[CH2+1]({carboxylic_carbon_idx}){hydroxy_pattern}")):
                        return True, "Contains two carboxylic groups with a hydroxy group on an alpha carbon"
    return False, "Does not have a hydroxy group on an alpha carbon to a carboxylic acid group"