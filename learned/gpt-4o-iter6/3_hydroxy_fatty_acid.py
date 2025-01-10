"""
Classifies: CHEBI:59845 3-hydroxy fatty acid
"""
from rdkit import Chem

def is_3_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acid based on its SMILES string.
    A 3-hydroxy fatty acid has a hydroxy functional group at the beta- or 3-position
    from the carboxylic acid and is characterized by a long carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 3-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Locate the carboxylic acid group using a SMARTS pattern
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxy_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxy_matches:
        return False, "No carboxylic acid group found"

    for carboxy_match in carboxy_matches:
        carboxylic_carbon_idx = carboxy_match[0]
        carboxylic_carbon = mol.GetAtomWithIdx(carboxylic_carbon_idx)

        # Set to check for hydroxy group at 3-position
        for neighbor in carboxylic_carbon.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:
                # First carbon (alpha position)
                alpha_carbon = neighbor
                for beta_neighbor in alpha_carbon.GetNeighbors():
                    if beta_neighbor.GetAtomicNum() == 6 and beta_neighbor.GetIdx() != carboxylic_carbon.GetIdx():
                        # Second carbon (beta position)
                        beta_carbon = beta_neighbor
                        for gamma_neighbor in beta_carbon.GetNeighbors():
                            if gamma_neighbor.GetAtomicNum() == 8:
                                # Check that the oxygen (hydroxy) is bonded to the gamma carbon
                                if any(h_neighbor.GetAtomicNum() == 1 for h_neighbor in gamma_neighbor.GetNeighbors()):
                                    return True, "Contains a hydroxy group at the 3-position"

    return False, "No hydroxy group at the 3-position found"