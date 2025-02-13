"""
Classifies: CHEBI:27208 unsaturated fatty acid
"""
from rdkit import Chem

def is_unsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acid based on its SMILES string.
    An unsaturated fatty acid contains at least one C=C or C#C bond, a carboxyl group, 
    and a sufficiently long aliphatic chain.
    
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

    # Check for carboxyl group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxyl group found (-COOH)"
    
    # Check for presence of unsaturation: C=C or C#C bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    triple_bond_pattern = Chem.MolFromSmarts("C#C")
    if not (mol.HasSubstructMatch(double_bond_pattern) or mol.HasSubstructMatch(triple_bond_pattern)):
        return False, "No unsaturation (C=C or C#C) found in the structure"

    # Check for long aliphatic chain; adjust chain length to a more suitable threshold
    carbon_chain_pattern = Chem.MolFromSmarts("[CH2]~[CH2]~[CH2]~[CH2]~[CH2]~[CH2]~[CH2]~[CH2]")  # At least 8 continuous carbons
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "Aliphatic chain is too short to be considered a fatty acid"

    return True, "Contains at least one unsaturation and a sufficiently long carbon chain with a carboxyl group, typical of unsaturated fatty acids"