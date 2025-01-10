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
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # General 3-hydroxy pattern accommodating stereochemistry and branches
    # (Carboxylic connection)-C-C-[C@H](O) or branching C
    three_hydroxy_pattern = Chem.MolFromSmarts("C(=O)O[C;R0][C;R0][C@H](O)")

    # Check matches for both simple and complex patterns
    if mol.HasSubstructMatch(three_hydroxy_pattern):
        return True, "Molecule matches 3-hydroxy fatty acid pattern"

    return False, "No 3-hydroxy group next to carboxyl on third carbon"

# Example use
example_smiles = "CCCCCCCCCC[C@@H](O)CC(O)=O"  # Example SMILES
result, reason = is_3_hydroxy_fatty_acid(example_smiles)
print(result, reason)