"""
Classifies: CHEBI:59845 3-hydroxy fatty acid
"""
from rdkit import Chem

def is_3_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acid based on its SMILES string.
    A 3-hydroxy fatty acid has a hydroxyl group on the 3rd carbon from
    the carboxyl carbon in its long aliphatic chain.

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

    # Pattern for a carboxylic acid connected to a carbon chain: C(=O)O
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Pattern to identify 3-hydroxy group in fatty acids
    # The pattern is represented as CC(CO)..C(=O)O where the third carbon has an OH
    three_hydroxy_fatty_acid_pattern = Chem.MolFromSmarts("CC(O)C(=O)O")
    if not mol.HasSubstructMatch(three_hydroxy_fatty_acid_pattern):
        return False, "No 3-hydroxy group next to carboxyl on third carbon"

    return True, "Molecule matches pattern for 3-hydroxy fatty acid"

# Example use
example_smiles = "CCCCCCCCCC[C@@H](O)CC(O)=O"  # Example SMILES
result, reason = is_3_hydroxy_fatty_acid(example_smiles)
print(result, reason)