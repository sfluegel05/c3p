"""
Classifies: CHEBI:27208 unsaturated fatty acid
"""
"""
Classifies: CHEBI:25681 unsaturated fatty acid
"""
from rdkit import Chem

def is_unsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acid based on its SMILES string.
    An unsaturated fatty acid is defined as any fatty acid containing at least one C=C or C#C bond.

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

    # Check for carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for at least one double (C=C) or triple (C#C) bond
    double_bond_pattern = Chem.MolFromSmarts("[CX3]=[CX3]")
    triple_bond_pattern = Chem.MolFromSmarts("[CX2]#[CX2]")
    
    has_double_bond = mol.HasSubstructMatch(double_bond_pattern)
    has_triple_bond = mol.HasSubstructMatch(triple_bond_pattern)
    
    if not (has_double_bond or has_triple_bond):
        return False, "No double or triple bonds found"

    return True, "Contains a carboxylic acid group and at least one double or triple bond"