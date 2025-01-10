"""
Classifies: CHEBI:26607 saturated fatty acid
"""
"""
Classifies: CHEBI:26607 saturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a saturated fatty acid based on its SMILES string.
    A saturated fatty acid is a carboxylic acid with a long, fully saturated carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a saturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if carboxylic_acid_pattern is None:
        return False, "Invalid SMARTS pattern for carboxylic acid"
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for long carbon chain (at least 6 carbons in the main chain)
    carbon_chain_pattern = Chem.MolFromSmarts("[CX4H2]~[CX4H2]~[CX4H2]~[CX4H2]~[CX4H2]~[CX4H2]")
    if carbon_chain_pattern is None:
        return False, "Invalid SMARTS pattern for carbon chain"
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "Carbon chain too short (less than 6 carbons)"

    # Check for unsaturation (double or triple bonds)
    unsaturated_pattern = Chem.MolFromSmarts("[CX3]=[CX3]|[CX2]#[CX2]")
    if unsaturated_pattern is None:
        return False, "Invalid SMARTS pattern for unsaturation"
    if mol.HasSubstructMatch(unsaturated_pattern):
        return False, "Carbon chain contains double or triple bonds (unsaturated)"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 6:
        return False, "Too few carbons for a fatty acid"
    if o_count < 2:
        return False, "Must have at least 2 oxygens (carboxylic acid group)"

    return True, "Contains a long, fully saturated carbon chain with a carboxylic acid group"