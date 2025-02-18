"""
Classifies: CHEBI:50760 N-hydroxy-alpha-amino-acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_hydroxy_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-hydroxy-alpha-amino-acid based on its SMILES string.
    An N-hydroxy-alpha-amino-acid has a carboxylic acid group and an N-hydroxy group on the alpha carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-hydroxy-alpha-amino-acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for carboxylic acid group
    carboxylic_acid = Chem.MolFromSmarts('[CX3](=[OX1])[OX2H1]')
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group"

    # Find all carboxylic acid carbons and check adjacent alpha carbons
    matches = mol.GetSubstructMatches(carboxylic_acid)
    for match in matches:
        carb_carbon_idx = match[0]
        carb_carbon = mol.GetAtomWithIdx(carb_carbon_idx)
        # Iterate through neighbors of carboxylic acid carbon (potential alpha carbons)
        for neighbor in carb_carbon.GetNeighbors():
            if neighbor.GetSymbol() == 'C':
                alpha_carbon = neighbor
                # Check neighbors of alpha carbon for nitrogen with hydroxyl group
                for alpha_neighbor in alpha_carbon.GetNeighbors():
                    if alpha_neighbor.GetSymbol() == 'N':
                        # Check if nitrogen has any hydroxyl groups attached
                        for n_neighbor in alpha_neighbor.GetNeighbors():
                            if n_neighbor.GetSymbol() == 'O' and n_neighbor.GetTotalNumHs() >= 1:
                                return True, "N-hydroxy group found on alpha carbon adjacent to carboxylic acid"

    return False, "No N-hydroxy group on alpha carbon adjacent to carboxylic acid"