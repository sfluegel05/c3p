"""
Classifies: CHEBI:27208 unsaturated fatty acid
"""
from rdkit import Chem

def is_unsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acid based on its SMILES string.
    An unsaturated fatty acid contains at least one C=C or C#C bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an unsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of carboxyl group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No terminal carboxyl group found (-COOH)"
    
    # Check for presence of unsaturation: C=C or C#C bonds
    unsaturation_pattern = Chem.MolFromSmarts("C=C | C#C")
    if not mol.HasSubstructMatch(unsaturation_pattern):
        return False, "No unsaturation (C=C or C#C) found in the structure"

    # Check for long aliphatic chain characteristic of fatty acids
    # Note: This is a simplistic check; a more sophisticated one would track the chain length.
    carbon_chain_pattern = Chem.MolFromSmarts("CCCC")
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "No long aliphatic carbon chain characteristic of fatty acids"

    return True, "Structure contains at least one unsaturation and terminal carboxyl group, typical of unsaturated fatty acids"