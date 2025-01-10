"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
"""
Classifies: CHEBI:17855 2-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxy fatty acid based on its SMILES string.
    A 2-hydroxy fatty acid has a hydroxy functional group in the alpha- or 2-position relative to the carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define patterns for carboxylic acid and 2-hydroxy group in the alpha position
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    hydroxy_alpha_pattern = Chem.MolFromSmarts("[C@H1,[C@@H1]]([O;H1])C(=O)O")  # including chiral center consideration

    # Check for carboxylic acid group
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for 2-hydroxy group in alpha position
    if not mol.HasSubstructMatch(hydroxy_alpha_pattern):
        return False, "No 2-hydroxy group found in the alpha position"
    
    # Ensure a long carbon chain (minimum length typical for fatty acids)
    alpah_carbon = Chem.MolFromSmarts("CCC(=O)O")
    if not mol.HasSubstructMatch(alpah_carbon):
        return False, "Insufficient length or wrong configuration in what should be alpha-chain for a fatty acid"
    
    # Count number of carbons beyond the mentioned minimum of fatty acids
    n_fatty_chain_carbons = rdMolDescriptors.CalcNumAliphaticCarbocycles(mol)
    if n_fatty_chain_carbons < 12:
        return False, "Insufficient carbon chain length for a fatty acid"

    return True, "Contains 2-hydroxy group in the alpha position of a fatty acid chain"