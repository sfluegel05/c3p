"""
Classifies: CHEBI:59845 3-hydroxy fatty acid
"""
from rdkit import Chem

def is_3_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acid based on its SMILES string.
    A 3-hydroxy fatty acid has a hydroxyl group on the 3rd carbon from
    the carboxyl carbon in its aliphatic chain.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for a carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxy_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxy_matches:
        return False, "No carboxylic acid group found"

    # Iterate over possible carboxylic acid positions
    for match in carboxy_matches:
        carboxylic_carbon_idx = match[0]  # Get the carboxylic carbon

        # Check for 3-hydroxy group at the third carbon from the carboxyl
        three_hydroxy_pattern = Chem.MolFromSmarts("C(=O)O[*][*]C(O)")
        if mol.HasSubstructMatch(three_hydroxy_pattern):
            return True, "Molecule matches 3-hydroxy fatty acid pattern"

    return False, "No 3-hydroxy group next to carboxyl on third carbon"

# Example use
example_smiles = "CCCCCCCCCCC[C@@H](O)CC(O)=O"  # Example SMILES
result, reason = is_3_hydroxy_fatty_acid(example_smiles)
print(result, reason)