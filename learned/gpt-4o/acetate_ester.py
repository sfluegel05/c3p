"""
Classifies: CHEBI:47622 acetate ester
"""
from rdkit import Chem

def is_acetate_ester(smiles: str):
    """
    Determines if a molecule is an acetate ester, specifically a carboxylic ester where
    the carboxylic acid component is acetic acid, based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an acetate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (False, "Invalid SMILES string")

    # SMARTS pattern for the acetate ester group
    # Acetate group: [C](=O)O[C], where the second [C] is a methyl group
    acetate_pattern = Chem.MolFromSmarts("C(=O)OC")

    # Look for acetate group
    if mol.HasSubstructMatch(acetate_pattern):
        matches = mol.GetSubstructMatches(acetate_pattern)
        for match in matches:
            # Verify that the ester oxygen is connected to a methyl group
            ester_oxygen = match[2]  # This is the oxygen in C(=O)O[C]
            methyl_carbon = match[3] # This should be part of CH3

            # Confirm that the carbon (methyl_carbon) has 3 hydrogens bonded to it, making it CH3
            carbon = mol.GetAtomWithIdx(methyl_carbon)
            if carbon.GetSymbol() == 'C' and sum(1 for nbr in carbon.GetNeighbors() if nbr.GetSymbol() == 'H') == 3:
                return (True, "Contains acetate ester group as part of acetic acid")

    return (False, "No acetate ester group found")